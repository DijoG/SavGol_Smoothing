##
require(tidyverse);require(gtable);require(see);require(grid);require(gridExtra)
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

##>------------------------ OA -------------------------------------
csvc <- list.files("D:/Savitzky/adatok/oa", pattern = "sok", full.names = T)
G <- read_csv2(csvc[1]) %>%
  mutate(Forest = "G")
M <- read_csv2(csvc[2]) %>%
  mutate(Forest = "M")
S <- read_csv2(csvc[3]) %>%
  mutate(Forest = "S")

OSAV <- bind_rows(list(G,M,S)) %>%
  mutate(Forest = factor(Forest),
         val = as.numeric(val),
         sav = ifelse(sav == "original", "b_0", sav)) %>%
  mutate(sav = factor(str_replace(sav, "^([^_]+_[^_]+)_(.*)$", "\\1::\\2")))

# Reverse the levels of sav
OSAV$sav <- factor(OSAV$sav, levels = rev(gtools::mixedsort(levels(OSAV$sav))))

# Calculate medians for each Forest and sav combination
medians <- OSAV %>%
  filter(sav == "b_0") %>%
  group_by(Forest) %>%
  summarise(median_val = median(val, na.rm = TRUE))

##> Visualisation
P <- OSAV %>%
  ggplot(aes(x = val, y = sav, fill = Forest)) +
  geom_boxplot(alpha = .5, outlier.size = 0.5, outlier.alpha = 0.3,
               position = position_dodge(width = 0.9), colour = "grey45",
               linewidth = .1) +
  geom_vline(data = medians, 
             aes(xintercept = median_val), color = "grey45", linewidth = .1, alpha = 1) +
  facet_wrap(~ Forest, ncol = 3) +
  labs(x = "Overall Accuracy", y = "Order::Window size") +
  scale_fill_manual(values = c("firebrick2", "forestgreen", "cyan2"),
                    name = "Forests") +
  scale_color_manual(values = c("firebrick2", "forestgreen", "cyan2"),
                     name = "Forests") +
  scale_x_continuous(breaks = c(0.7, 0.8, 0.9), labels = c("0.7", "0.8", "0.9"),
                     expand = expansion(mult = c(0, 0))) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.key.size = unit(5, "mm"),
        strip.text = element_text(size = 9, margin = margin(2, 0, 2, 0))) +
  guides(fill = guide_legend(override.aes = list(alpha = .6), ncol = 5))

P

##>---- Polynomial order vs. window size (forest sites) and OA -----
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

# Generate/wrangle data for plotting
OSAV_WO <- 
  OSAV %>%
  mutate(Window = as.numeric(str_extract(sav, "(?<=b_)\\d+")),
         Order = as.numeric(str_extract(sav, "(?<=::)\\d+"))) %>%
  mutate(Order = ifelse(is.na(Order), 0, Order))

OSAV_WOP <- 
  OSAV_WO %>%
  group_by(Forest) %>%
  mutate(
    baseline_val = mean(val[sav == "b_0"], na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Order, Forest) %>%
  mutate(accuracy_diff = (val - baseline_val) * 100) %>%
  ungroup() %>%
  filter(sav != "b_0") 

OSAV_plot <- 
  OSAV_WOP %>%
  mutate(OrderForest = str_c(Order, "_", Forest)) %>%
  group_by(OrderForest, Window, Forest, Order) %>%
  summarise(mean_accuracy_diff = mean(accuracy_diff)) %>%
  ungroup()

# Factorize Order and Window
OSAV_plot$Order <- factor(OSAV_plot$Order, levels = sort(unique(OSAV_plot$Order)))
OSAV_plot$Window <- factor(OSAV_plot$Window, levels = sort(unique(OSAV_plot$Window)))

##> Plot
ggplot(OSAV_plot, aes(x = Window, y = Order, fill = mean_accuracy_diff)) +
  geom_tile(color = NA) +
  geom_text(aes(label = round(mean_accuracy_diff, 1)), 
            size = 2.5, 
            color = "grey30") +
  facet_wrap(~ Forest) +  # Important: allows different y-axis labels per facet
  scale_fill_gradient2(low = "red",
                       mid = "white", 
                       high = "darkgreen",
                       midpoint = 0,
                       name = "Accuracy difference from baseline (%)",
                       guide = guide_colorbar(title.position = "right",  
                                              title.hjust = 0.5,         
                                              barwidth = unit(0.3, "cm"),
                                              barheight = unit(6, "cm"))) +
  scale_x_discrete(breaks = unique(OSAV_plot$Window), expand = expansion(mult = c(0.0, 0.0))) +
  scale_y_discrete(breaks = seq(1, 17, 1), expand = expansion(mult = c(0.0, 0.0))) +
  labs(x = "Window size",y = "Order") +
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 15, margin = margin(2, 0, 2, 0)),
        axis.text = element_text(size = 9),
        axis.title = element_text(),
        legend.position = "right",
        legend.justification = "left",  
        legend.margin = margin(t = 0, r = 0, b = 0, l = -5),  
        legend.title = element_text(size = 10, colour = "grey25", angle = 90),
        legend.text = element_text(size = 8))

##>---------------------- VRR analysis -----------------------------
FF <- read_csv("D:/Savitzky/adatok/oa/all_sites.csv") %>%
  mutate(sav = ifelse(sav == "original", "b_0", sav)) %>%
  mutate(sav = factor(str_replace(sav, "^([^_]+_[^_]+)_(.*)$", "\\1::\\2")),
         site = case_when(
           site == "Mecsek" ~ "M",
           site == "Gemenc" ~ "G",
           TRUE ~ "S"
         )) %>%
  rename(Forest = site) 
FF$sav <- factor(FF$sav, levels = rev(gtools::mixedsort(levels(FF$sav))))

highlight_points <- data.frame(
  Forest = c("G", "M", "S"),
  sav = c("b_5::1", "b_3::1", "b_5::1"),
  label = c("b_5::1", "b_3::1", "b_5::1")
)

highlight_df <- FF %>%
  inner_join(highlight_points, by = c("Forest", "sav"))

##> Visualisation
FF %>%
  ggplot(aes(x = mean_vrr, y = val)) +
  geom_point(aes(size = std_vrr, color = W_minus_P), alpha = 0.5) +
  scale_color_gradient2(low = "dodgerblue3", 
                        mid = "forestgreen", 
                        high = "yellow",
                        midpoint = 10) +
  scale_size_continuous(breaks = c(0.05, 0.10, 0.15),
                        range = c(1, 5),
                        labels = c("< 0.05", "0.05 - 0.1", "> 0.1")) +    
  facet_wrap(~Forest, scales = "free_x") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.1)))+
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.04)))+
  geom_label(data = highlight_df,
             aes(x = mean_vrr, y = val, label = label),
             hjust = -0.05, vjust = -0.3, size = 4,
             label.size = 0,           # Remove border
             label.padding = unit(0.2, "lines"),  # Minimal padding
             fill = "grey100",
             color = "firebrick2",
             alpha = 0) +
  # Optional: add a point outline to make them stand out
  geom_point(data = highlight_df,
             aes(x = mean_vrr, y = val),
             shape = 21, size = 2, stroke = 2, color = "firebrick2", fill = NA) +
  labs(x = "Mean Variance Reduction Ratio",
       y = "Median Overall Accuracy",
       size = "Smoothing inconsistency (Std VRR)",
       color = "Polynomial fitting imbalance (W-P)") +
  guides(size = guide_legend(direction = "horizontal", title.position = "top",
                             override.aes = list(color = "grey90")),
         color = guide_colorbar(barwidth = 8, barheight = .7, alpha = .9)) + 
  theme(legend.position = "top",           
        legend.box = "horizontal",
        legend.key.size = unit(5, "mm"),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 13),
        strip.text = element_text(size = 12, margin = margin(2, 0, 2, 0)),
        legend.spacing.x = unit(1.3, "cm"),
        legend.box.spacing = unit(0.1, "cm"))  

##>---------------------- Smoothing Imagery ------------------------
require(terra);require(tidyterra)

tifs <- list.files("D:/Savitzky/adatok/gemenc_fig5", pattern = ".tif$", full.names = T)
r <- rast(tifs[3])
r2 <- rast(tifs[2])
r3 <- rast(tifs[1])

# Extract wavelengths from band names
band_names <- names(r)

# Function to extract wavelength from the parentheses at the end
extract_wavelength <- function(name) {
  # Extract the number inside the last parentheses
  matches <- regmatches(name, gregexpr("\\([0-9.]+[^\\)]*\\)", name))
  last_match <- matches[[1]][length(matches[[1]])]
  # Remove parentheses and convert to numeric
  as.numeric(gsub("\\(|\\)", "", last_match))
}

wavelengths <- sapply(band_names, extract_wavelength)

# Create band_info data frame
band_info <- data.frame(
  band_index = 1:length(band_names),
  band_name = band_names,
  wavelength_nm = wavelengths,
  min_value = minmax(r)[1,],
  max_value = minmax(r)[2,],
  range = minmax(r)[2,] - minmax(r)[1,]
)

# Sort by wavelength to see the full spectrum
band_info_sorted <- band_info[order(band_info$wavelength_nm), ]
print("All bands sorted by wavelength:")
print(band_info_sorted[, c("band_index", "wavelength_nm")])

# Find RGB bands
# Red: 620-750 nm
red_candidates <- band_info[band_info$wavelength_nm >= 620 & band_info$wavelength_nm <= 750, ]
# Green: 495-570 nm
green_candidates <- band_info[band_info$wavelength_nm >= 495 & band_info$wavelength_nm <= 570, ]
# Blue: 450-495 nm
blue_candidates <- band_info[band_info$wavelength_nm >= 450 & band_info$wavelength_nm <= 495, ]

# Display candidates
cat("\n🔴 RED BANDS (620-750 nm):\n")
if(nrow(red_candidates) > 0) {
  print(red_candidates[, c("band_index", "wavelength_nm")])
  # Select the closest to 670 nm (typical red)
  red_idx <- which.min(abs(red_candidates$wavelength_nm - 670))
  red_band <- red_candidates$band_index[red_idx]
  red_wl <- red_candidates$wavelength_nm[red_idx]
  cat("→ Selected Red Band:", red_band, sprintf("(%.3f nm)\n", red_wl))
} else {
  cat("No red bands found\n")
  red_band <- NULL
}

# 59

cat("\n🟢 GREEN BANDS (495-570 nm):\n")
if(nrow(green_candidates) > 0) {
  print(green_candidates[, c("band_index", "wavelength_nm")])
  # Select the closest to 550 nm (typical green)
  green_idx <- which.min(abs(green_candidates$wavelength_nm - 550))
  green_band <- green_candidates$band_index[green_idx]
  green_wl <- green_candidates$wavelength_nm[green_idx]
  cat("→ Selected Green Band:", green_band, sprintf("(%.3f nm)\n", green_wl))
} else {
  cat("No green bands found\n")
  green_band <- NULL
}

# 24

cat("\n🔵 BLUE BANDS (450-495 nm):\n")
if(nrow(blue_candidates) > 0) {
  print(blue_candidates[, c("band_index", "wavelength_nm")])
  # Select the closest to 470 nm (typical blue)
  blue_idx <- which.min(abs(blue_candidates$wavelength_nm - 470))
  blue_band <- blue_candidates$band_index[blue_idx]
  blue_wl <- blue_candidates$wavelength_nm[blue_idx]
  cat("→ Selected Blue Band:", blue_band, sprintf("(%.3f nm)\n", blue_wl))
} else {
  cat("No blue bands found\n")
  blue_band <- NULL
}

# 1

##> Plot

# Load your data and reproject
tifs <- list.files("D:/Savitzky/adatok/gemenc_fig5", pattern = ".tif$", full.names = T)
r <- rast(tifs[3]) %>%
  terra::project("EPSG:23700")  # gemenc_o.tif
r2 <- rast(tifs[2]) %>%
  terra::project("EPSG:23700")  # gemenc_5.tif  
r3 <- rast(tifs[1]) %>%
  terra::project("EPSG:23700")  # gemenc_19.tif

# Set your RGB bands
red_band <- 59
green_band <- 24
blue_band <- 1

# Create RGB composites
rgb1 <- c(r[[red_band]], r[[green_band]], r[[blue_band]])
rgb2 <- c(r2[[red_band]], r2[[green_band]], r2[[blue_band]])
rgb3 <- c(r3[[red_band]], r3[[green_band]], r3[[blue_band]])

# PDF for vector quality
pdf("RGB_composites_abc.pdf", width = 15, height = 5)

par(mfrow = c(1, 3),
    mar = c(3, 3, 2, 1),
    oma = c(0, 0, 0, 0))

# Plot (a) with fixed coordinates for EPSG:23700
plotRGB(rgb1, r = 1, g = 2, b = 3, stretch = "lin", axes = F)
rect(639486, 96002, 639525, 96030, col = "white", border = NA)  
text(x = 639502, y = 96017, labels = "(a)", cex = 3.3, font = 1, col = "grey25")

# Plot (b)
plotRGB(rgb2, r = 1, g = 2, b = 3, stretch = "lin", axes = F)
rect(639486, 96002, 639525, 96030, col = "white", border = NA)
text(x = 639502, y = 96017, labels = "(b)", cex = 3.3, font = 1, col = "grey25")

# Plot (c)
plotRGB(rgb3, r = 1, g = 2, b = 3, stretch = "lin", axes = F)
rect(639486, 96002, 639525, 96030, col = "white", border = NA)
text(x = 639502, y = 96017, labels = "(c)", cex = 3.3, font = 1, col = "grey25")

dev.off()
