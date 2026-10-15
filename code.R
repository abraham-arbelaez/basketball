# Space Jammed - CHANCE figures and analysis
# Data: shots-2019.csv, players.csv
# X is across the court, Y points toward half court (feet).
# Rim at (0, 0); court spans X = +/-25 and Y = -5.25 to 41.75.
# Figures go to figures/, model fits are cached in cache/.

suppressPackageStartupMessages({
  library(tidyverse)
  library(spatstat)
  library(mgcv)
  library(sf)
  library(spdep)
  library(CARBayes)
  library(patchwork)
})

dir.create("figures", showWarnings = FALSE)
dir.create("cache", showWarnings = FALSE)
set.seed(2019)

res_file <- "figures/results.txt"
cat("Space Jammed -- numbers quoted in the text\n", file = res_file)
report <- function(...) {
  line <- paste0(...)
  cat(line, "\n")
  cat(line, "\n", file = res_file, append = TRUE)
}

# 1. Data
raw <- read.csv("shots-2019.csv")[, -1]

# Parse the play descriptions instead of splitting on spaces.
# Names with suffixes or multiple words caused problems before.
m <- str_match(raw$play, "<br>(.+?) (made|missed) ([23])-pointer from (\\d+) ft")
stopifnot(!anyNA(m[, 1]))

shots <- raw %>%
  transmute(
    game_id, team, quarter, time_remaining,
    player  = m[, 2],
    made    = as.integer(m[, 3] == "made"),
    value   = as.integer(m[, 4]),             # 2 or 3
    dist_pbp = as.integer(m[, 5]),            # distance reported in play-by-play
    xraw    = as.numeric(sub("px", "", x)) / 10,  # ft from baseline (marker corner)
    yraw    = as.numeric(sub("px", "", y)) / 10   # ft across the court
  )

# Estimate rim location from reported shot distances (rounded to feet).
ok <- shots$dist_pbp < 30
rim <- optim(c(5, 25), function(p)
  sum((sqrt((shots$xraw[ok] - p[1])^2 + (shots$yraw[ok] - p[2])^2) - shots$dist_pbp[ok])^2))$par

shots <- shots %>%
  mutate(X = yraw - rim[2],
         Y = xraw - rim[1],
         r = sqrt(X^2 + Y^2),
         theta = atan2(Y, abs(X)) * 180 / pi)   # 0 = baseline/corner, 90 = straight on

# Three-point boundary
beyond_arc <- function(X, Y) ifelse(Y < sqrt(23.75^2 - 22^2), abs(X) > 22, sqrt(X^2 + Y^2) > 23.75)
agree <- mean(beyond_arc(shots$X, shots$Y) == (shots$value == 3))

# Leave out backcourt shots
n_all <- nrow(shots)
shots <- filter(shots, Y <= 41.75, abs(X) <= 25, Y >= -5.25)

# Some names differ between the shot and roster files
name_fix <- c("J.J. Redick" = "JJ Redick", "Danuel House" = "Danuel House Jr.",
              "Sviatoslav Mykhailiuk" = "Svi Mykhailiuk", "James Ennis" = "James Ennis III",
              "Otto Porter" = "Otto Porter Jr.", "Marko Guduric" = "Marko Gudurić",
              "Juan Hernangómez" = "Juancho Hernangómez", "Enes Kanter" = "Enes Freedom",
              "Charlie Brown" = "Charlie Brown Jr.", "Nicolas Claxton" = "Nic Claxton",
              "Lonnie Walker" = "Lonnie Walker IV", "Wesley Iwundu" = "Wes Iwundu")
players <- read_csv("players.csv", show_col_types = FALSE) %>% dplyr::select(Player, Pos)
shots <- shots %>%
  mutate(player_std = coalesce(name_fix[player], player)) %>%
  left_join(players, by = c("player_std" = "Player")) %>%
  mutate(Pos = factor(Pos, levels = c("PG", "SG", "SF", "PF", "C")))
stopifnot(!anyNA(shots$Pos))

n_games <- n_distinct(shots$game_id)
team_games <- 2 * n_games

report("Shots in file: ", n_all, "; in half-court window: ", nrow(shots),
       " (", n_all - nrow(shots), " back-court heaves dropped)")
report("Games: ", n_games, " (", team_games, " team-games); players: ", n_distinct(shots$player))
report(sprintf("Rim location in raw basketball-reference coords: (%.2f, %.2f) ft", rim[1], rim[2]))
report(sprintf("Official 3-pt line reproduces the play-by-play 2/3 label for %.2f%% of shots", 100 * agree))
report(sprintf("League FG%%: %.1f%%; points per shot (excl. FTs): %.3f",
               100 * mean(shots$made), mean(shots$made * shots$value)))
report(sprintf("Share of attempts that are 3s: %.1f%%", 100 * mean(shots$value == 3)))
report("FG% by position: ",
       paste(sprintf("%s %.1f%%", levels(shots$Pos), 100 * tapply(shots$made, shots$Pos, mean)), collapse = ", "))

# Court lines
arc <- function(cx, cy, r, from, to, n = 200) {
  t <- seq(from, to, length.out = n)
  tibble(X = cx + r * cos(t), Y = cy + r * sin(t))
}
corner_y <- sqrt(23.75^2 - 22^2)          # where the corner 3 meets the arc
court_lines <- bind_rows(
  tibble(X = c(-25, 25, 25, -25, -25), Y = c(-5.25, -5.25, 41.75, 41.75, -5.25)) %>% mutate(id = "boundary"),
  tibble(X = c(-8, -8, 8, 8), Y = c(-5.25, 13.75, 13.75, -5.25)) %>% mutate(id = "paint"),
  arc(0, 13.75, 6, 0, pi) %>% mutate(id = "ft_top"),
  arc(0, 0, 0.75, 0, 2 * pi) %>% mutate(id = "rim"),
  tibble(X = c(-3, 3), Y = c(-1.25, -1.25)) %>% mutate(id = "board"),
  bind_rows(tibble(X = -4, Y = -1.25), arc(0, 0, 4, pi, 0), tibble(X = 4, Y = -1.25)) %>% mutate(id = "restricted"),
  bind_rows(tibble(X = -22, Y = -5.25),
            arc(0, 0, 23.75, pi - atan2(corner_y, 22), atan2(corner_y, 22)),
            tibble(X = 22, Y = -5.25)) %>% mutate(id = "three"),
  arc(0, 41.75, 6, pi, 2 * pi) %>% mutate(id = "center")
)
ft_bottom <- arc(0, 13.75, 6, pi, 2 * pi) %>% mutate(id = "ft_bottom")

court <- function(colour = "grey25", lw = 0.3, ylim = c(-5.25, 41.75), three = TRUE) {
  cl <- if (three) court_lines else filter(court_lines, id != "three")
  list(
    geom_path(data = cl, aes(X, Y, group = id), colour = colour, linewidth = lw, inherit.aes = FALSE),
    geom_path(data = ft_bottom, aes(X, Y), colour = colour, linewidth = lw, linetype = "22", inherit.aes = FALSE),
    coord_fixed(xlim = c(-25, 25), ylim = ylim, expand = FALSE, clip = "on")
  )
}

theme_court <- function(base = 9) {
  theme_void(base_size = base) +
    theme(plot.title = element_text(face = "bold", size = base + 1, hjust = 0),
          plot.subtitle = element_text(colour = "grey30", hjust = 0),
          plot.tag = element_text(face = "bold"),
          legend.position = "bottom",
          legend.title = element_text(size = base - 1, hjust = 0.5),
          legend.title.position = "top",
          legend.text = element_text(size = base - 2),
          legend.key.height = unit(0.32, "cm"),
          legend.key.width = unit(1.1, "cm"),
          plot.margin = margin(4, 4, 4, 4))
}
theme_set(theme_court())

save_fig <- function(p, name, w, h) {
  ggsave(file.path("figures", paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path("figures", paste0(name, ".png")), p, width = w, height = h, dpi = 400, bg = "white")
}

YLIM <- c(-5.25, 32)   # most figures only need the front court

# Shot chart
# Jitter for the plot only; otherwise the coordinate grid is very visible.
# Use the original locations for the models.
fig1 <- ggplot(shots, aes(X, Y)) +
  geom_jitter(size = 0.35, alpha = 0.07, stroke = 0, shape = 16, width = 0.25, height = 0.25) +
  court(colour = "#C8102E", lw = 0.35) +
  theme_court()
save_fig(fig1, "fig1_shot_chart", 4.5, 4.4)

# Point pattern / shot intensity
W <- owin(c(-25, 25), c(-5.25, 41.75))
pp <- ppp(shots$X, shots$Y, window = W, marks = shots$Pos)
pix <- 0.25                                # 3-inch pixels

# CV bandwidths were almost zero because the coordinates repeat so often.
# Using 0.75 ft here; keep Scott's bandwidth in the results for reference.
sigma <- 0.75
report(sprintf("Kernel bandwidth: sigma = %.2f ft (Scott's rule would give %.2f x %.2f ft)",
               sigma, bw.scott(pp)[1], bw.scott(pp)[2]))

lam <- density(unmark(pp), sigma = sigma, eps = pix, edge = TRUE, positive = TRUE)
lam_df <- as.data.frame(lam) %>% rename(X = x, Y = y) %>%
  mutate(per_game = value / team_games)       # shots per team-game per sq ft

# Sidebar: points, counts, smoothed intensity
top_player <- names(sort(table(shots$player), decreasing = TRUE))[1]
ps <- filter(shots, player == top_player)
pp1 <- ppp(ps$X, ps$Y, window = W)
q <- 2
cnt <- ps %>% mutate(Xc = floor(X / q) * q + q / 2, Yc = floor(Y / q) * q + q / 2) %>%
  count(Xc, Yc)
lam1 <- as.data.frame(density(pp1, sigma = 1.5, eps = pix, edge = TRUE, positive = TRUE)) %>%
  rename(X = x, Y = y)
report(sprintf("Sidebar player: %s, %d shots", top_player, nrow(ps)))

sb_a <- ggplot(ps, aes(X, Y)) +
  geom_point(size = 0.35, alpha = 0.5, colour = "grey10", stroke = 0) +
  court(ylim = YLIM) + labs(title = "1. Shots as points")
sb_b <- ggplot(cnt, aes(Xc, Yc, fill = n)) +
  geom_tile(width = q, height = q) +
  scale_fill_viridis_c(option = "magma", direction = -1, trans = "sqrt", name = "Shots per\n2x2 ft square") +
  court(ylim = YLIM) + labs(title = "2. Count them in a grid")
sb_c <- ggplot(lam1, aes(X, Y, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma", direction = -1, trans = "sqrt", name = "Shots per sq ft\n(kernel estimate)") +
  court(ylim = YLIM) + labs(title = "3. Shrink the squares & smooth")
sb_theme <- theme(legend.key.width = unit(0.6, "cm"), legend.title = element_text(size = 7))
sidebar1 <- (sb_a + sb_theme) | (sb_b + sb_theme) | (sb_c + sb_theme)
save_fig(sidebar1, "sidebar1_intensity", 7.2, 3.1)

# League intensity (log scale)
brks <- c(1e-4, 1e-3, 1e-2, 1e-1, 1)
fig2 <- ggplot(lam_df, aes(X, Y, fill = pmax(per_game, 1e-4))) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma", direction = -1, trans = "log10",
                       breaks = brks, labels = c("0.0001", "0.001", "0.01", "0.1", "1"),
                       limits = c(1e-4, max(lam_df$per_game)), oob = scales::squish,
                       name = "Shots per team-game\nper square foot (log scale)") +
  court(colour = "grey40", ylim = YLIM) +
  theme_court()
save_fig(fig2, "fig2_intensity", 4.2, 4.0)

# Shot locations by position, relative to the league
# rho_k(s) = [lambda_k(s) / n_k] / [lambda(s) / n] ; 1 = same shot profile as league
sig_rel <- 1.5
lam_all_s <- density(unmark(pp), sigma = sig_rel, eps = pix, edge = TRUE, positive = TRUE)
rel_df <- map_dfr(levels(shots$Pos), function(k) {
  ppk <- unmark(pp[marks(pp) == k])
  lk <- density(ppk, sigma = sig_rel, eps = pix, edge = TRUE, positive = TRUE)
  rr <- (lk / npoints(ppk)) / (lam_all_s / npoints(pp))
  as.data.frame(rr) %>% mutate(Pos = k, dens = as.data.frame(lam_all_s)$value)
}) %>%
  rename(X = x, Y = y) %>%
  mutate(Pos = factor(Pos, levels = levels(shots$Pos),
                      labels = c("Point guards", "Shooting guards", "Small forwards",
                                 "Power forwards", "Centers")),
         value = ifelse(dens / team_games < 2e-4, NA, value))   # mask near-empty floor

fig3 <- ggplot(rel_df, aes(X, Y, fill = log2(value))) +
  geom_raster() +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-2, 2), oob = scales::squish,
                       breaks = -2:2, labels = c("¼×", "½×", "same", "2×", "4×"),
                       na.value = "grey92",
                       name = "Shot intensity relative to league average") +
  court(colour = "grey30", lw = 0.25, ylim = YLIM) +
  facet_wrap(~Pos, nrow = 1) +
  theme_court() +
  theme(strip.text = element_text(face = "bold", size = 8, margin = margin(b = 3)),
        legend.key.width = unit(1.4, "cm"))
save_fig(fig3, "fig3_relative_intensity_position", 7.2, 2.25)

# ---------------------------------------------------------------------------
# Spatial GAM for P(make | X, Y)
# Try a few GP ranges and keep the fit with the lowest fREML score.
# ---------------------------------------------------------------------------
gam_file <- "cache/gam_gp.rds"
if (!file.exists(gam_file)) {
  ranges <- c(0.75, 1, 1.5, 2, 3, 5, 10)
  fits <- map(ranges, function(rho) {
    t0 <- Sys.time()
    f <- bam(made ~ s(X, Y, bs = "gp", k = 400, m = c(3, rho)),
             family = binomial, data = shots, method = "fREML", discrete = TRUE, nthreads = 4)
    message(sprintf("range %4.1f ft: fREML %.1f, AIC %.1f (%.0fs)", rho, f$gcv.ubre, AIC(f),
                    as.numeric(Sys.time() - t0, units = "secs")))
    f
  })
  prof <- tibble(range = ranges, fREML = map_dbl(fits, "gcv.ubre"), AIC = map_dbl(fits, AIC))
  best <- fits[[which.min(prof$fREML)]]
  saveRDS(list(fit = best, profile = prof), gam_file)
}
gg <- readRDS(gam_file)
gam_fit <- gg$fit
report("GP range profile (Matern 3/2, k = 400):")
for (i in seq_len(nrow(gg$profile)))
  report(sprintf("   range %4.1f ft: fREML %.1f  AIC %.1f", gg$profile$range[i], gg$profile$fREML[i], gg$profile$AIC[i]))
report(sprintf("Chosen GP range: %.0f ft; edf = %.1f", gg$profile$range[which.min(gg$profile$fREML)],
               sum(gam_fit$edf) - 1))

grid <- expand_grid(X = seq(-25, 25, by = 0.25), Y = seq(-5.25, 41.75, by = 0.25))
pr <- predict(gam_fit, newdata = grid, type = "link", se.fit = TRUE)
grid <- grid %>%
  mutate(eta = pr$fit, se = pr$se.fit,
         p = plogis(eta),
         lo = plogis(eta - 1.96 * se), hi = plogis(eta + 1.96 * se),
         value = ifelse(beyond_arc(X, Y), 3, 2),
         xp = p * value)
# Mask areas with very few attempts
dens_at <- interp.im(lam, grid$X, grid$Y)
grid$support <- dens_at / team_games > 2e-4

league_fg <- mean(shots$made)
league_pps <- mean(shots$made * shots$value)

g_show <- filter(grid, support, Y <= YLIM[2])
fig4a <- ggplot(g_show, aes(X, Y, fill = p)) +
  geom_raster() +
  scale_fill_viridis_c(option = "mako", direction = -1, labels = scales::percent,
                       limits = c(0.25, 0.70), oob = scales::squish,
                       name = "Estimated FG%") +
  court(ylim = YLIM) + labs(title = "Probability of a make")
fig4b <- ggplot(g_show, aes(X, Y, fill = 100 * (hi - lo))) +
  geom_raster() +
  scale_fill_viridis_c(option = "rocket", direction = -1, name = "95% interval width\n(percentage points)") +
  court(ylim = YLIM) + labs(title = "How sure are we?")
fig4c <- ggplot(g_show, aes(X, Y, fill = xp)) +
  geom_raster() +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(0.6, 1.5), oob = scales::squish,
                       values = scales::rescale(c(0.6, league_pps, 1.5)),
                       name = "Expected points\nper shot") +
  court(ylim = YLIM) + labs(title = "Expected points")
narrow_key <- theme(legend.key.width = unit(0.7, "cm"))
fig4 <- (fig4a + narrow_key) | (fig4b + narrow_key) | (fig4c + narrow_key)
save_fig(fig4, "fig4_continuous_surface", 7.2, 2.95)

# Residual variogram: aggregate to 2-ft cells and adjust for distance.
# Standardize by estimated binomial variance before computing pairwise differences.
cell <- 2
cells <- shots %>%
  mutate(Xc = floor(X / cell) * cell + cell / 2, Yc = floor(Y / cell) * cell + cell / 2,
         rc = sqrt(Xc^2 + Yc^2)) %>%
  group_by(Xc, Yc, rc) %>% summarise(n = n(), makes = sum(made), .groups = "drop") %>%
  filter(n >= 20, rc >= 4, rc <= 30)
trend <- gam(cbind(makes, n - makes) ~ s(rc, k = 20), family = binomial, data = cells)
p_tr <- fitted(trend)
cells$z <- (cells$makes / cells$n - p_tr) / sqrt(p_tr * (1 - p_tr) / cells$n)
dmat <- as.matrix(dist(cells[, c("Xc", "Yc")]))
hv <- outer(cells$z, cells$z, "-")^2 / 2
ut <- upper.tri(dmat)
vg <- tibble(h = dmat[ut], g = hv[ut]) %>%
  filter(h <= 20) %>%
  mutate(bin = cut(h, breaks = seq(1, 21, by = 2))) %>%
  filter(!is.na(bin)) %>%
  group_by(bin) %>% summarise(h = mean(h), gamma = mean(g), npairs = n(), .groups = "drop")
report(sprintf("Variogram (2-ft cells, n = %d): gamma %.2f at %.1f ft rising to %.2f beyond 10 ft",
               nrow(cells), vg$gamma[1], vg$h[1], mean(vg$gamma[vg$h > 10])))
sidebar2 <- ggplot(vg, aes(h, gamma)) +
  geom_hline(yintercept = 1, linetype = "22", colour = "grey50") +
  annotate("text", x = 20, y = 1, label = "pure chance", vjust = -0.5, hjust = 1, size = 2.8, colour = "grey40") +
  geom_line(colour = "grey60") + geom_point(size = 2, colour = "#1F4E79") +
  scale_x_continuous(limits = c(0, 20.5), expand = c(0, 0)) +
  expand_limits(y = 0.9) +
  labs(x = "Distance between two 2-ft squares (ft)",
       y = "Semivariance",
       title = "Nearby spots shoot alike",
       subtitle = "FG% residuals after removing distance to rim") +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(colour = "grey30"),
        panel.grid.minor = element_blank())
save_fig(sidebar2, "sidebar2_variogram", 3.6, 2.8)

# Hexagonal lattice / Leroux CAR
court_poly <- st_polygon(list(rbind(c(-25, -5.25), c(25, -5.25), c(25, 32), c(-25, 32), c(-25, -5.25)))) %>%
  st_sfc()
hex <- st_make_grid(court_poly, cellsize = 2, square = FALSE) %>%
  st_sf(geometry = .) %>%
  filter(st_within(st_centroid(geometry), court_poly, sparse = FALSE)[, 1]) %>%
  mutate(hex_id = row_number())
shots_sf <- st_as_sf(shots, coords = c("X", "Y"), remove = FALSE)
shots_sf$hex_id <- st_intersects(shots_sf, hex) %>% map_int(~ if (length(.x)) .x[1] else NA_integer_)

hex_dat <- shots_sf %>% st_drop_geometry() %>% filter(!is.na(hex_id)) %>%
  group_by(hex_id) %>%
  summarise(n = n(), makes = sum(made), share3 = mean(value == 3), .groups = "drop")
hex <- hex %>% inner_join(hex_dat, by = "hex_id") %>%
  mutate(raw_fg = makes / n, value = 2 + share3)

nb <- poly2nb(hex, snap = 1e-6)
keep <- card(nb) > 0
hex <- hex[keep, ]
nb <- poly2nb(hex, snap = 1e-6)
Wmat <- nb2mat(nb, style = "B")

mor <- moran.test(hex$raw_fg, nb2listw(nb, style = "W"))
report(sprintf("Lattice: %d hexagons (2-ft), median %d shots/hex (range %d-%d)",
               nrow(hex), median(hex$n), min(hex$n), max(hex$n)))
report(sprintf("Moran's I of raw hex FG%%: %.3f (p = %.2g)", mor$estimate[1], mor$p.value))

car_file <- "cache/car.rds"
if (!file.exists(car_file)) {
  car <- S.CARleroux(makes ~ 1, trials = hex$n, family = "binomial", data = st_drop_geometry(hex),
                     W = Wmat, burnin = 20000, n.sample = 120000, thin = 20, verbose = FALSE)
  saveRDS(car, car_file)
}
car <- readRDS(car_file)
rho_s <- car$samples$rho
report(sprintf("Leroux CAR: posterior mean rho = %.3f (95%% CrI %.3f-%.3f); tau2 = %.3f",
               mean(rho_s), quantile(rho_s, .025), quantile(rho_s, .975), mean(car$samples$tau2)))
report("CAR fit summary:"); capture.output(print(car$summary.results)) %>% walk(report)

hex <- hex %>% mutate(car_fg = car$fitted.values / n,
                      car_pps = car_fg * value,
                      shrink = car_fg - raw_fg)

lim_fg <- c(0.25, 0.70)
hex_theme <- list(court(ylim = YLIM), theme_court())
fig5a <- ggplot(hex) + geom_sf(aes(fill = raw_fg), colour = NA) + hex_theme +
  scale_fill_viridis_c(option = "mako", direction = -1, labels = scales::percent, limits = lim_fg,
                       oob = scales::squish, name = "FG%") +
  coord_sf(xlim = c(-25, 25), ylim = YLIM, expand = FALSE) +
  labs(title = "Raw FG% in each hexagon")
fig5b <- ggplot(hex) + geom_sf(aes(fill = car_fg), colour = NA) + hex_theme +
  scale_fill_viridis_c(option = "mako", direction = -1, labels = scales::percent, limits = lim_fg,
                       oob = scales::squish, name = "FG%") +
  coord_sf(xlim = c(-25, 25), ylim = YLIM, expand = FALSE) +
  labs(title = "Borrowing from neighbors (CAR)")
fig5c <- ggplot(hex) + geom_sf(aes(fill = car_pps), colour = NA) + hex_theme +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(0.6, 1.5), oob = scales::squish,
                       values = scales::rescale(c(0.6, league_pps, 1.5)),
                       name = "Expected points\nper shot") +
  coord_sf(xlim = c(-25, 25), ylim = YLIM, expand = FALSE) +
  labs(title = "Expected points per shot")
fig5 <- (fig5a + narrow_key) | (fig5b + narrow_key) | (fig5c + narrow_key)
suppressMessages(save_fig(fig5, "fig5_lattice_car", 7.2, 2.95))

# Check how much CAR estimates differ from raw proportions
fig5d <- ggplot(st_drop_geometry(hex), aes(n, 100 * shrink)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_point(alpha = 0.6, size = 1.2, colour = "#1F4E79") +
  scale_x_log10() +
  labs(x = "Shots in hexagon (log scale)", y = "CAR estimate minus raw FG%\n(percentage points)",
       title = "Small samples borrow the most") +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
save_fig(fig5d, "fig5d_shrinkage", 3.6, 2.8)

# ---------------------------------------------------------------------------
# Recovering the three-point line from shooting locations
# The known line is used for comparison, not for estimating its location.
# ---------------------------------------------------------------------------

# Gradient of log intensity
lam_fine <- density(unmark(pp), sigma = 0.75, eps = 0.2, edge = TRUE, positive = TRUE)
ll <- log(lam_fine)
M <- as.matrix(ll)
dX <- (cbind(M[, -1], NA) - cbind(NA, M[, -ncol(M)])) / (2 * ll$xstep)
dY <- (rbind(M[-1, ], NA) - rbind(NA, M[-nrow(M), ])) / (2 * ll$ystep)
grad <- ll; grad$v <- sqrt(dX^2 + dY^2)
grad_df <- tibble(X = as.vector(rasterx.im(grad)), Y = as.vector(rastery.im(grad)),
                  value = as.vector(grad$v),
                  dens = as.vector(lam_fine$v) / team_games) %>%
  mutate(value = ifelse(dens < 5e-4 | Y > 30, NA, value))

fig6a <- ggplot(grad_df, aes(X, Y, fill = value)) +
  geom_raster() +
  scale_fill_gradient(low = "white", high = "black", limits = c(0, 2), oob = scales::squish,
                      na.value = "white", guide = "none") +
  geom_path(data = filter(court_lines, id == "boundary"), aes(X, Y), inherit.aes = FALSE,
            colour = "grey50", linewidth = 0.3) +
  coord_fixed(xlim = c(-25, 25), ylim = YLIM, expand = FALSE) +
  labs(title = "Where shot density changes sharply") +
  theme_court()

# For each angle, look for the largest increase in log intensity
lam_line <- density(unmark(pp), sigma = 0.5, eps = 0.1, edge = TRUE, positive = TRUE)
r_grid <- seq(15, 30, by = 0.1)
angles <- seq(1.5, 88.5, by = 3)
est_line <- map_dfr(angles, function(a) {
  map_dfr(c(-1, 1), function(side) {
    Xs <- side * r_grid * cos(a * pi / 180); Ys <- r_grid * sin(a * pi / 180)
    inside <- Ys <= 41.75 & abs(Xs) <= 25
    v <- rep(NA_real_, length(r_grid))
    v[inside] <- log(interp.im(lam_line, Xs[inside], Ys[inside]))
    dv <- c(NA, diff(v)) / 0.1
    j <- which.max(dv)
    tibble(angle = a, side = side, r_hat = r_grid[j], X = Xs[j], Y = Ys[j])
  })
})
true_r <- function(a) ifelse(23.75 * sin(a * pi / 180) < corner_y, 22 / cos(a * pi / 180), 23.75)
est_line <- mutate(est_line, r_true = true_r(angle), err = r_hat - r_true)
report(sprintf("Data-driven 3-pt line: mean abs. error %.2f ft (max %.2f ft) over %d directions",
               mean(abs(est_line$err)), max(abs(est_line$err)), nrow(est_line)))

fig6b <- ggplot() +
  court(colour = "grey70", lw = 0.3, ylim = YLIM) +
  geom_point(data = est_line, aes(X, Y), colour = "#C8102E", size = 1.3) +
  labs(title = "The three-point line, rediscovered") +
  theme_court()

# Radial density, FG%, and points per shot
zone <- function(theta, Y) ifelse(Y < corner_y, "Corner", "Above the break")
prof <- shots %>%
  filter(r >= 10, r <= 30) %>%
  mutate(zone = zone(theta, Y), rb = floor(r * 2) / 2 + 0.25) %>%
  group_by(zone, rb) %>%
  summarise(n = n(), fg = mean(made), pps = mean(made * value), .groups = "drop")
# Approximate area of each radial bin on a fine grid
ag <- expand_grid(X = seq(-24.95, 24.95, by = 0.1), Y = seq(-5.2, 41.7, by = 0.1)) %>%
  mutate(r = sqrt(X^2 + Y^2), theta = atan2(Y, abs(X)) * 180 / pi, zone = zone(theta, Y),
         rb = floor(r * 2) / 2 + 0.25) %>%
  count(zone, rb, name = "cells") %>% mutate(area = cells * 0.01)
prof <- left_join(prof, ag, by = c("zone", "rb")) %>%
  mutate(dens = n / area / team_games) %>%
  complete(zone, rb) %>%
  # Skip FG% and PPS where the bin is too small
  mutate(thin = is.na(n) | n < 200,
         fg = ifelse(thin, NA, fg), pps = ifelse(thin, NA, pps))

zcols <- c("Corner" = "#C8102E", "Above the break" = "#1F4E79")
vlines <- tibble(zone = c("Corner", "Above the break"), r = c(22, 23.75))
prof_theme <- list(
  geom_vline(data = vlines, aes(xintercept = r, colour = zone), linetype = "22", linewidth = 0.4),
  geom_line(linewidth = 0.6), geom_point(size = 0.6),
  scale_colour_manual(values = zcols, name = NULL),
  theme_minimal(base_size = 9),
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 9))
)
fig6c1 <- ggplot(prof, aes(rb, dens, colour = zone)) + prof_theme +
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) +
  labs(x = NULL, y = NULL, title = "Shots per team-game per sq ft (log)")
fig6c2 <- ggplot(prof, aes(rb, fg, colour = zone)) + prof_theme +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Distance from rim (ft)", y = NULL, title = "Field-goal percentage")
fig6c3 <- ggplot(prof, aes(rb, pps, colour = zone)) + prof_theme +
  labs(x = NULL, y = NULL, title = "Points per shot")

fig6 <- (fig6a | fig6b) + plot_annotation(tag_levels = "a")
save_fig(fig6, "fig6_discover_basketball", 6.2, 2.85)

fig7 <- (fig6c1 | fig6c2 | fig6c3) + plot_layout(guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom"))
save_fig(fig7, "fig7_radial_profiles", 7.2, 2.7)

report("Radial profile just inside vs just outside the line (team-game density, FG%):")
prof %>%
  filter((zone == "Corner" & rb %in% c(21.25, 21.75, 22.25, 22.75)) |
         (zone == "Above the break" & rb %in% c(22.75, 23.25, 24.25, 24.75))) %>%
  transmute(line = sprintf("   %-15s r=%5.2f  dens=%.4f  FG=%.1f%%  PPS=%.2f  (n=%d)",
                           zone, rb, dens, 100 * fg, pps, n)) %>%
  pull(line) %>% walk(report)

report("\nsessionInfo:"); capture.output(sessionInfo()) %>% head(12) %>% walk(report)
