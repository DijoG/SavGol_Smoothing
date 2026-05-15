##
require(INLA);require(tidyverse);require(gtable);require(see);require(grid);require(gridExtra)
##

##> Theme
theme_set(theme_lucid(base_family = "Lato") + # theme_lucid
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.line.y = element_blank(),
                  axis.line.x = element_blank(), # "grey98"
                  panel.background = element_rect(fill = "grey98"), # color = NA
                  axis.ticks = element_blank(),
                  axis.text.x = element_text(size = 10),
                  axis.text.y = element_text(size = 10),
                  axis.title = element_text(size = 14),
                  strip.text = element_text(size = 8, margin = margin(2, 0, 2, 0)),
                  strip.background = element_rect(color = "transparent")))

##> Load Data and wrangle a bit
OSAV <- read_csv2("D:/Savitzky/adatok/SMG_fafaj_oa.csv") %>%
  mutate(sav = factor(sav),
         Species = Species_e,
         Species = factor(Species),
         Forest = factor(Forest),
         Forest_Species = factor(str_c(Forest, "_", Species)),         # interaction 1
         sav_Forest_Species = factor(str_c(sav, "_", Forest_Species)), # interaction 2
         sav_Forest = factor(str_c(sav, "_", Forest)),                 # interaction 3
         val = as.numeric(val)) %>%
  mutate(sav = factor(sav, levels = gtools::mixedsort(levels(sav))))

# Parameters from the POV of the original 
OSAVo <- OSAV %>% filter(sav == "b_0")
orig_MEAN <- OSAVo$val %>% mean()
orig_PRECI <- 1/(sd(OSAVo$val)^2)

# Formula
FORMULA <- val ~ 1 +
  f(sav, model = "iid") +
  f(Forest, model = "iid") +
  f(Species, model = "iid") +
  f(Forest_Species, model = "iid") +                                  # interaction 1
  f(sav_Forest_Species, model = "iid") +                              # interaction 2
  f(sav_Forest, model = "iid")                                        # interaction 3

##----------------
##> Run INLA magic
##----------------
# names(inla.models()$likelihood) -> no student-t, so let's stick to a modified (heavy-tailed) Gaussian
SAVinla <- 
  inla(
    FORMULA, 
    family = "gaussian", 
    data = OSAV, 
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE, 
                           waic = TRUE, 
                           config = TRUE,
                           return.marginals.predictor=TRUE),
    control.fixed = list(mean.intercept = orig_MEAN, 
                         prec.intercept = orig_PRECI),
    # Weaker prior for precision (1, 0.02) -> (0.5, 0.005) for heavier tails on the likelihood ->
    control.family = list(hyper = list(prec = list(prior = "loggamma", param = c(0.5, 0.005)))),  
    control.inla = list(strategy = "adaptive"))

summary(SAVinla)                           # very low DIC and WAIC indicate a good model fit

## -------------------------------------------------------------------------------
## -------------------- FIXED EFFECTS
# Extract posterior predictive metrics 
preds <- 
  OSAV %>%
  bind_cols(SAVinla$summary.fitted.values %>%
              as.data.frame()) %>%
  mutate(RESIDUALS = val - mean) %>%
  rename(Predicted = "mean",
         Observed = "val")

preds %>%
  filter(sav == "b_0") %>%
  pull(Observed) %>% mean() == orig_MEAN
pred_MEAN <- 
  preds %>%
  filter(sav == "b_0") %>%
  pull(Predicted) %>% mean() 
pred_MEAN - orig_MEAN                       # difference is really small;)

## -------------------------------------------------------------------------------
## -------------------- POSTERIOR PREDICTIVE CHECKS (model-fitting check)
# 01 
preds %>%
  ggplot(aes(x = RESIDUALS)) +
  geom_density(fill = "grey35", color = NA) + 
  geom_jitter(aes(y = -3), width = 0.1, height = 1, alpha = 0.02, color = "grey25", shape = 16) +
  geom_vline(xintercept = 0, col = "firebrick2", linewidth = .2) +
  labs(subtitle = "Posterior Predictive Check") +
  labs(x = "Residuals") +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.01))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme(plot.subtitle = element_text(hjust = 0.5, colour = "grey45"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 9),
        axis.title.y = element_blank())

# 02
preds %>%
  ggplot(aes(x = Predicted, y = Observed)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "firebrick2") +
  labs(subtitle = "Posterior Predictive Check",
       x = "Predicted OA",
       y = "Observed OA") 

flab <- 
  preds %>%
  group_by(Forest) %>%
  summarise(flab = unique(Forest), .groups = "drop") %>%
  mutate(Observed = 0.03)

# 03
preds %>%
  ggplot(aes(x = Observed, y = Forest)) +
  geom_point(aes(color = "Observed", alpha = .1), shape = 16) +  
  geom_point(aes(x = Predicted, y = Forest, color = "Predicted", alpha = .05), shape = 16) +  
  facet_wrap(~ sav, ncol = 7) +
  labs(x = "Overall Accuracy", color = "Type",
       subtitle = "Posterior Predictive Check") +
  geom_text(data = flab, aes(x = Observed, y = Forest, label = flab), 
            inherit.aes = FALSE, hjust = -0.1, size = 3, col = "grey65") +
  scale_y_discrete(expand = expansion(mult = c(0.1, 0.1))) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), 
                     expand = expansion(mult = c(0.0, 0.0))) +
  scale_color_manual(values = c("Observed" = "firebrick2", "Predicted" = "forestgreen")) +  
  scale_alpha_identity() + 
  theme(legend.position = c(.85, .03),
        legend.background = element_rect(fill = "grey98", colour = "transparent"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 7),
        plot.subtitle = element_text(hjust = 0.5, colour = "grey45"),
        legend.key.size = unit(6, "mm"),
        legend.title = element_blank()) +
  #legend.box.spacing = margin(2)) +
  guides(colour = guide_legend(override.aes = list(alpha = .7)))


## ---------------------------------------------------------------------------------
## -------------------- RANDOM EFFECTS - Updated version with 95% credible intervals
# Extract random effects from INLA output
# 1) Summary
get_RANDOMstatEFF <- function(inlamodel) {
  params = names(inlamodel$summary.random)
  all_eff = 
    map2_dfr(inlamodel$summary.random, params, ~ 
               mutate(.x, effect = .y))
  
  stat_EFF = 
    all_eff %>%
    filter(0 > `0.975quant` | 0 < `0.025quant`)
  
  non_stat_EFF = 
    anti_join(all_eff, stat_EFF, by = colnames(all_eff)) %>%
    arrange(abs(mean))
  
  return(list(stat_EFF = stat_EFF,
              non_stat_EFF = non_stat_EFF))
}
random_effects_sum <- get_RANDOMstatEFF(SAVinla)
random_effects_sum$stat_EFF
random_effects_sum$non_stat_EFF 


# 2) Marginals
get_RANDOMmarginals <- function(inlamodel, data) {
  
  params = names(SAVinla$marginals.random)
  
  OUT = list()
  for (i in seq_along(params)) {
    effm = SAVinla$marginals.random[[params[i]]]
    dat = OSAV[[params[i]]] %>% levels()
    
    OUT[[i]] = 
      map2_df(effm, dat, ~ data.frame(inla.smarginal(.x)) %>%
                mutate(effect = .y,
                       group = params[i]))
  }
  out = bind_rows(OUT) 
  return(out)
}
random_effects <- get_RANDOMmarginals(SAVinla, OSAV) %>%
  mutate(effect = factor(effect)) %>%
  mutate(effect = factor(effect, levels = gtools::mixedsort(levels(effect))))

## non stat EFF sav = no difference between the performance of orig and smoothed
nstatEFF <- str_extract(random_effects_sum$non_stat_EFF$ID, "^b_\\d+::\\d+$") 
nstateffect <- data.frame(effect = nstatEFF[!is.na(nstatEFF)],
                          x = 0,
                          y = 2,
                          group = "sav")
nstatEFF <- data.frame(effect = random_effects %>%
                         filter(group == "sav") %>%
                         pull(effect) %>% unique() %>%
                         as.character(),
                       x = NA,
                       y = NA,
                       group = "sav")
nstatEFF <- nstatEFF %>%
  left_join(nstateffect, by = "effect", suffix = c("", ".new")) %>%
  mutate(x = ifelse(is.na(x.new), x, x.new),
         y = ifelse(is.na(y.new), y, y.new)) %>%
  select(effect, x, y, group)

# Extract marginals and compute credible intervals
get_RANDOMmarginals_withCI <- function(inlamodel, data) {
  
  params = names(inlamodel$marginals.random)
  
  OUT = list()
  for (i in seq_along(params)) {
    effm = inlamodel$marginals.random[[params[i]]]
    dat = data[[params[i]]] %>% levels()
    
    # Calculate 95% credible intervals for each effect
    ci_list <- list()
    for (j in seq_along(effm)) {
      marginal <- effm[[j]]
      # Compute credible interval
      ci <- inla.qmarginal(c(0.025, 0.975), marginal)
      ci_list[[j]] <- data.frame(
        lower = ci[1],
        upper = ci[2],
        median = inla.qmarginal(0.5, marginal),
        mean = inla.emarginal(function(x) x, marginal),
        effect = dat[j],
        group = params[i]
      )
    }
    
    OUT[[i]] <- bind_rows(ci_list)
  }
  return(bind_rows(OUT))
}

# Get credible intervals for random effects
random_effects_ci <- get_RANDOMmarginals_withCI(SAVinla, OSAV) %>%
  mutate(effect = factor(effect)) %>%
  mutate(effect = factor(effect, levels = gtools::mixedsort(levels(effect))))

## non stat EFF sav = no difference between the performance of orig and smoothed
nstatEFF <- str_extract(random_effects_sum$non_stat_EFF$ID, "^b_\\d+::\\d+$") 
nstateffect <- data.frame(effect = nstatEFF[!is.na(nstatEFF)],
                          x = 0,
                          y = 2,
                          group = "sav")
nstatEFF <- data.frame(effect = random_effects_ci %>%
                         filter(group == "sav") %>%
                         pull(effect) %>% unique() %>%
                         as.character(),
                       x = NA,
                       y = NA,
                       group = "sav")
nstatEFF <- nstatEFF %>%
  left_join(nstateffect, by = "effect", suffix = c("", ".new")) %>%
  mutate(x = ifelse(is.na(x.new), x, x.new),
         y = ifelse(is.na(y.new), y, y.new)) %>%
  select(effect, x, y, group)

## 01: Plot random effects "sav" with 95% credible interval bars
p_ci_bars <- 
  random_effects_ci %>%
  filter(group == "sav") %>%
  ggplot(aes(x = mean, y = effect)) +
  geom_vline(xintercept = 0, col = "firebrick2", linewidth = 0.1, alpha = 0.5) +
  geom_point(size = 1, color = "grey35") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), 
                 height = 0.2, linewidth = 0.3, color = "grey35") +
  geom_text(nstatEFF, mapping = aes(x, y, label = "*"), 
            col = "firebrick2", size = 5, vjust = 0.8) +
  labs(x = "Posterior effect (mean with 95% CI)", y = "Smoothing parameter") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  theme_minimal(base_family = "Lato") +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.title = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey98", color = NA),
    axis.line.x = element_line(color = "grey80")
  )

# 02: Density plots with shaded credible intervals
# First get the full marginals for plotting
random_effects_full <- get_RANDOMmarginals(SAVinla, OSAV) %>%
  mutate(effect = factor(effect)) %>%
  mutate(effect = factor(effect, levels = gtools::mixedsort(levels(effect))))

# Join with CI data
random_effects_with_ci <- random_effects_full %>%
  left_join(random_effects_ci, by = c("effect", "group"))

# Plot densities with shaded 95% credible shaded intervals
p_ci_shaded <- 
  random_effects_with_ci %>%
  filter(group == "sav") %>%
  ggplot(aes(x = x, y = y)) +
  geom_vline(aes(xintercept = 0, col = "0.0"), linewidth = .1) +
  # Shade the 95% credible interval
  geom_rect(data = . %>% distinct(effect, .keep_all = TRUE),
            aes(xmin = lower, xmax = upper, ymin = 0, ymax = Inf),
            fill = "grey85", alpha = 0.3, inherit.aes = FALSE) +
  # Add density line
  geom_line(col = "grey35", linewidth = 0.3) +
  # Add mean point
  geom_point(data = . %>% distinct(effect, .keep_all = TRUE),
             aes(x = mean, y = 0.5), 
             size = 1, color = "firebrick2", shape = 16) +
  geom_text(nstatEFF, mapping = aes(x, y, label = "*"), 
            col = "firebrick2", size = 5) +
  labs(x = "Posterior effect") +
  facet_wrap(~ effect, ncol = 7) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 11),
    legend.position = c(.86, .03),
    legend.background = element_rect(fill = "grey98", colour = "transparent"), 
    legend.title = element_blank(),
    legend.key.width = unit(1, "mm"),
    panel.background = element_rect(fill = "grey98", color = NA),
    strip.text = element_text(size = 7),
    strip.background = element_rect(fill = "grey90", color = NA),
    panel.spacing = unit(0.2, "lines")
  )

# 03: Compact summary plot (all effects in one panel with different colors)
p_compact <- 
  random_effects_ci %>%
  filter(group == "sav") %>%
  mutate(effect_num = as.numeric(str_extract(effect, "\\d+")),
         is_significant = ifelse(lower > 0 | upper < 0, "Significant", "Non-significant")) %>%
  ggplot(aes(x = mean, y = reorder(effect, effect_num), color = is_significant)) +
  geom_vline(xintercept = 0, col = "grey60", linewidth = 0.2, linetype = "dashed") +
  geom_point(size = 1.5) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, linewidth = 0.5) +
  scale_color_manual(values = c("Significant" = "firebrick2", "Non-significant" = "grey35")) +
  labs(x = "Posterior effect (mean with 95% CI)", 
       y = "Smoothing parameter",
       color = "95% CI excludes 0") +
  theme_minimal(base_family = "Lato") +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    axis.text = element_text(size = 8)
  )

# Display the plots
p_ci_bars    # Simple interval bars
p_ci_shaded  # Densities with shaded intervals
p_compact    # Compact summary

# ONLY SHADED ONE IS NEEDED!
g <- ggplotGrob(p_ci_shaded)
strip_positions <- which(g$layout$name %in% grep("strip-t", g$layout$name, value = TRUE))

if (length(strip_positions) > 0) {
  first_strip_index = strip_positions[1]
  first_strip = g$grobs[[first_strip_index]]
  print(first_strip)
  strip_label = "*"  
} else {
  stop("No facet strips found!")
}

# PLOT >
grid.newpage()
grid.draw(g)  

grid.text(
  label = strip_label,  
  x = unit(0.0878, "npc"),  
  y = unit(0.98, "npc"),  
  gp = gpar(col = "grey45", cex = 1, fontface = "bold", family = "Lato black")
)

## ----------------------------------------------------------------------------------
## ------------------- stat EFF "Forests"::"Species" = influencing OA "significantly"
# Get credible intervals for Forest_Species
forest_species_ci <- random_effects_ci %>%
  filter(group == "Forest_Species") %>%
  mutate(effect = str_replace(effect, "_", "::"))

statEFF <- random_effects_sum$stat_EFF %>%
  filter(effect == "Forest_Species")
stateffect <- data.frame(effect = statEFF[!is.na(statEFF)],
                         x = random_effects_ci %>%
                           filter(effect == "b_0") %>%
                           pull(mean) %>% mean,
                         y = 2,
                         group = "Forest_Species")
statEFF <- data.frame(effect = random_effects_ci %>%
                        filter(group == "Forest_Species") %>%
                        pull(effect) %>% unique() %>%
                        as.character(),
                      x = NA,
                      y = NA,
                      group = "Forest_Species")
statEFF <- statEFF %>%
  left_join(stateffect, by = "effect", suffix = c("", ".new")) %>%
  mutate(x = ifelse(is.na(x.new), x, x.new),
         y = ifelse(is.na(y.new), y, y.new)) %>%
  select(effect, x, y, group) %>%
  mutate(effect = str_replace(effect, "_", "::"))

## Plot random effects "Forest" :: "Species" interaction with shaded 95% CI
p_forest_species <- 
  random_effects_with_ci %>%
  filter(group == "Forest_Species") %>%
  mutate(effect = str_replace(effect, "_", "::")) %>%
  ggplot(aes(x, y)) +
  # Vertical line at 0
  geom_vline(aes(xintercept = 0, col = "0.0"), linewidth = .1) +
  # Shade the 95% credible interval
  geom_rect(data = . %>% distinct(effect, .keep_all = TRUE),
            aes(xmin = lower, xmax = upper, ymin = 0, ymax = Inf),
            fill = "grey85", alpha = 0.3, inherit.aes = FALSE) +
  # Add density line
  geom_line(col = "grey35", linewidth = .1) +
  # Add mean point
  geom_point(data = . %>% distinct(effect, .keep_all = TRUE),
             aes(x = mean, y = 0.2), 
             size = 1, color = "firebrick2", shape = 16) +
  # Add asterisk for significant effects
  geom_text(statEFF, mapping = aes(x, y, label = "*"), col = "firebrick2", size = 5) +
  scale_color_manual(values = "firebrick2", labels = "0.0") +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  scale_x_continuous(expand = expansion(mult = c(0.0, 0.0))) +
  labs(x = "Posterior effect") +
  facet_wrap(~ effect, ncol = 4) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 11),
    legend.position = c(.76, .1),
    legend.title = element_blank(),
    legend.key.width = unit(1, "mm"),
    legend.background = element_rect(fill = "grey98", colour = "transparent"),
    strip.text = element_text(size = 9),
    panel.background = element_rect(fill = "grey98", color = NA)
  )

## Two-panel plot: Forest and Species with shaded 95% CI
# Get credible intervals for Forest and Species
forest_ci <- random_effects_ci %>% filter(group == "Forest")
species_ci <- random_effects_ci %>% filter(group == "Species")

# Get full marginals for Forest and Species
forest_marginals <- random_effects_with_ci %>% filter(group == "Forest")
species_marginals <- random_effects_with_ci %>% filter(group == "Species")

# Get b_0 mean for reference
b0_mean <- random_effects_ci %>%
  filter(effect == "b_0") %>%
  pull(mean) %>% mean()

# Panel A: Forest effects
A <- 
  forest_marginals %>%
  ggplot(aes(x, y)) +
  # Vertical line at b_0 mean
  geom_vline(xintercept = b0_mean, col = "firebrick2", linewidth = .1) +
  # Shade the 95% credible interval
  geom_rect(data = forest_ci,
            aes(xmin = lower, xmax = upper, ymin = 0, ymax = Inf),
            fill = "grey85", alpha = 0.3, inherit.aes = FALSE) +
  # Add density line
  geom_line(col = "grey35", linewidth = .1) +
  # Add mean point
  geom_point(data = forest_ci,
             aes(x = mean, y = 0.1), 
             size = 1, color = "firebrick2", shape = 16) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  scale_x_continuous(expand = expansion(mult = c(0.0, 0.0))) +
  labs(x = "Posterior effect", tag = "(a)") +
  facet_wrap(~ effect, ncol = 3) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 11),
    plot.tag = element_text(colour = "grey25", size = 13, hjust = .5),
    plot.tag.position = "topright",
    strip.text = element_text(size = 9),
    panel.background = element_rect(fill = "grey98", color = NA)
  )

# Panel B: Species effects
B <- 
  species_marginals %>%
  ggplot(aes(x, y)) +
  # Vertical line at b_0 mean
  geom_vline(aes(xintercept = b0_mean, col = "0.0"), linewidth = .1) +
  # Shade the 95% credible interval
  geom_rect(data = species_ci,
            aes(xmin = lower, xmax = upper, ymin = 0, ymax = Inf),
            fill = "grey85", alpha = 0.3, inherit.aes = FALSE) +
  # Add density line
  geom_line(col = "grey35", linewidth = .1) +
  # Add mean point
  geom_point(data = species_ci,
             aes(x = mean, y = 0.3), 
             size = 1, color = "firebrick2", shape = 16) +
  scale_color_manual(values = "firebrick2", labels = "0.0") +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  scale_x_continuous(expand = expansion(mult = c(0.0, 0.0))) +
  labs(x = "Posterior effect", tag = "(b)") +
  facet_wrap(~ effect, ncol = 3) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 11),
    legend.position = c(.83, .07),
    legend.background = element_rect(fill = "grey98", colour = "transparent"), 
    legend.title = element_blank(),
    legend.key.width = unit(1, "mm"),
    plot.tag = element_text(colour = "grey25", size = 13, hjust = .5),
    plot.tag.position = "topright",
    strip.text = element_text(size = 9),
    panel.background = element_rect(fill = "grey98", color = NA)
  )

# Arrange the two panels
grid.arrange(A, B, ncol = 1, heights = c(1, 4))

# Forest_Species interaction 
p_forest_species

## ----------------------------------------------------------------------------------
## ------------------------------------- Best model obtained from Marginals (summary)
SAVM <- 
  random_effects_sum$non_stat_EFF %>%
  filter(effect == "sav") %>%
  mutate(MEAN = pred_MEAN + mean,
         Q_0025 = pred_MEAN + `0.025quant`,
         Q_0975 = pred_MEAN + `0.975quant`) %>%
  arrange(MEAN) %>%
  filter(!ID == "b_0") 

SAVM <- 
  SAVM %>%
  mutate(ID = factor(ID, levels = SAVM$ID))

modosit <- function(labels, start = 1) {
  labels[seq(start, length(labels), by = 2)] = ""
  return(labels)
}

calc_element("axis.text.y", theme_get()) # -> "grey50"

p1 <-
  ggplot(SAVM, aes(y = ID, x = MEAN, group = 1)) +
  geom_vline(xintercept = round(pred_MEAN, 2), size = .1, col = "grey65") +
  geom_ribbon(aes(xmin = Q_0025, xmax = Q_0975), fill = "grey75", alpha = 0.2) +  
  geom_line(color = "firebrick2", size = .2) + 
  geom_point(color = "firebrick2", size = 3, shape = 16, alpha = .5) +
  labs(y = "Order::Window size", x = "Overall Accuracy") +
  scale_x_continuous(expand = expansion(mult = c(0.0, 0.0)), 
                     breaks = c(0.78, 0.8, 0.83, 0.86, 0.88),
                     labels = c("0.78", "0.80", "0.83", "0.86","0.88"),
                     limits = c(0.767, 0.89)) +
  scale_y_discrete(expand = expansion(mult = c(0.00, 0.00)),
                   labels = modosit(levels(SAVM$ID), start = 2)) +
  theme(axis.text.y = element_text(angle = 0, size = 10, margin = margin(r = 5, unit = "pt")),
        axis.text.x = element_text(size = 9),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        panel.spacing = unit(0, "lines"))
p2 <-
  ggplot(SAVM, aes(y = ID, x = MEAN, group = 1)) +
  scale_y_discrete(
    expand = expansion(mult = c(0.0, 0.0)),
    labels = modosit(levels(SAVM$ID), start = 1)) +
  theme(
    axis.text.y = element_text(angle = 0, size = 10, color = "grey50", hjust = 0, vjust = .4, margin = margin(r = -30, unit = "pt")),
    axis.text.x = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.spacing = unit(0, "lines"))

require(patchwork)
layout <- "
AAAAAAAAB"
p1 + p2 + plot_layout(design = layout)


