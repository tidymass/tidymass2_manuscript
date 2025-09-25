library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(ggplot2)
library(tidyverse)

load("3_data_analysis/6_supplementary_note/object_in_house_annotation.rda")
load("2_data/ms2_original/mpsnyder_rplc_ms2.rda")
load("3_data_analysis/6_supplementary_note/fpa_validation.rda")

urine_annotation <- object_new@annotation_table

ms2_data <- object_new@ms2_data[[1]]

# for M363T580_POS
which(ms2_data@variable_id == "M363T580_POS")

spec_1 <- ms2_data@ms2_spectra[[3160]] %>% as.data.frame()

spec_1_standard <- mpsnyder_rplc_ms2@spectra.data[["Spectra.positive"]][["RPLC_679"]][["NCE25"]]


# normalization
spec_1 <- spec_1 %>%
  mutate(rel_intensity = intensity / max(intensity, na.rm = TRUE))
spec_1_standard <- spec_1_standard %>%
  mutate(rel_intensity = intensity / max(intensity, na.rm = TRUE))

# add labels
spec_1$source <- "Sample"
spec_1_standard$source <- "Standard"
spec_1_standard$rel_intensity <- -spec_1_standard$rel_intensity

# combine
spec_all <- bind_rows(spec_1, spec_1_standard)

## calculate axis limits
mz_min <- floor(min(spec_all$mz, na.rm = TRUE))
mz_max <- ceiling(max(spec_all$mz, na.rm = TRUE))
x_pad  <- (mz_max - mz_min) * 0.02
x_lim  <- c(mz_min - x_pad, mz_max + x_pad)

p <- ggplot(spec_all, aes(x = mz, y = rel_intensity, color = source)) +  
  
  geom_segment(aes(xend = mz, yend = 0), linewidth = 0.53) +  # 0.35 * 1.5 ≈ 0.53
  scale_color_manual(values = c("Standard" = "red", "Sample" = "black")) +
  guides(color = "none") +
  labs(
    x = "Mass to charge ratio (m/z)",
    y = "Relative intensity"
  ) +
  
  scale_y_continuous(
    limits = c(-1.1, 1.1),  
    breaks = seq(-1, 1, 0.5),
    labels = function(y) abs(y)
  ) +
  coord_cartesian(xlim = x_lim, expand = FALSE) +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.ticks.length = unit(2.5, "pt"),
    plot.title = element_text(hjust = 0, face = "bold"),
    plot.subtitle = element_text(hjust = 0),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.25),
    panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60")

print(p)

ggsave("4_manuscript/Supplementary_note/fpa_ms2/CE25/M363T580_POS.pdf", plot = p, width = 8, height = 6, dpi = 300)




# for M311T602_POS
which(ms2_data@variable_id == "M311T602_POS")

spec_1 <- ms2_data@ms2_spectra[[2603]] %>% as.data.frame()

spec_1_standard <- mpsnyder_rplc_ms2@spectra.data[["Spectra.positive"]][["RPLC_403"]][["NCE50"]]


# normalization
spec_1 <- spec_1 %>%
  mutate(rel_intensity = intensity / max(intensity, na.rm = TRUE))
spec_1_standard <- spec_1_standard %>%
  mutate(rel_intensity = intensity / max(intensity, na.rm = TRUE))

# add labels
spec_1$source <- "Sample"
spec_1_standard$source <- "Standard"
spec_1_standard$rel_intensity <- -spec_1_standard$rel_intensity

# combine
spec_all <- bind_rows(spec_1, spec_1_standard)

## calculate axis limits
mz_min <- floor(min(spec_all$mz, na.rm = TRUE))
mz_max <- ceiling(max(spec_all$mz, na.rm = TRUE))
x_pad  <- (mz_max - mz_min) * 0.02
x_lim  <- c(mz_min - x_pad, mz_max + x_pad)

p <- ggplot(spec_all, aes(x = mz, y = rel_intensity, color = source)) +  
  
  geom_segment(aes(xend = mz, yend = 0), linewidth = 0.53) +  # 0.35 * 1.5 ≈ 0.53
  scale_color_manual(values = c("Standard" = "red", "Sample" = "black")) +
  guides(color = "none") +
  labs(
    x = "Mass to charge ratio (m/z)",
    y = "Relative intensity"
  ) +
  
  scale_y_continuous(
    limits = c(-1.1, 1.1),  
    breaks = seq(-1, 1, 0.5),
    labels = function(y) abs(y)
  ) +
  coord_cartesian(xlim = x_lim, expand = FALSE) +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.ticks.length = unit(2.5, "pt"),
    plot.title = element_text(hjust = 0, face = "bold"),
    plot.subtitle = element_text(hjust = 0),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.25),
    panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60")

print(p)

ggsave("4_manuscript/Supplementary_note/fpa_ms2/CE50/M311T602_POS.pdf", plot = p, width = 8, height = 6, dpi = 300)


# for M274T390_POS
which(ms2_data@variable_id == "M274T390_POS")

spec_1 <- ms2_data@ms2_spectra[[2081]] %>% as.data.frame()

spec_1_standard <- mpsnyder_rplc_ms2@spectra.data[["Spectra.positive"]][["RPLC_259"]][["NCE25"]]


# normalization
spec_1 <- spec_1 %>%
  mutate(rel_intensity = intensity / max(intensity, na.rm = TRUE))
spec_1_standard <- spec_1_standard %>%
  mutate(rel_intensity = intensity / max(intensity, na.rm = TRUE))

# add labels
spec_1$source <- "Sample"
spec_1_standard$source <- "Standard"
spec_1_standard$rel_intensity <- -spec_1_standard$rel_intensity

# combine
spec_all <- bind_rows(spec_1, spec_1_standard)

## calculate axis limits
mz_min <- floor(min(spec_all$mz, na.rm = TRUE))
mz_max <- ceiling(max(spec_all$mz, na.rm = TRUE))
x_pad  <- (mz_max - mz_min) * 0.02
x_lim  <- c(mz_min - x_pad, mz_max + x_pad)

p <- ggplot(spec_all, aes(x = mz, y = rel_intensity, color = source)) +  
  
  geom_segment(aes(xend = mz, yend = 0), linewidth = 0.53) +  # 0.35 * 1.5 ≈ 0.53
  scale_color_manual(values = c("Standard" = "red", "Sample" = "black")) +
  guides(color = "none") +
  labs(
    x = "Mass to charge ratio (m/z)",
    y = "Relative intensity"
  ) +
  
  scale_y_continuous(
    limits = c(-1.1, 1.1),  
    breaks = seq(-1, 1, 0.5),
    labels = function(y) abs(y)
  ) +
  coord_cartesian(xlim = x_lim, expand = FALSE) +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.ticks.length = unit(2.5, "pt"),
    plot.title = element_text(hjust = 0, face = "bold"),
    plot.subtitle = element_text(hjust = 0),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.25),
    panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60")

print(p)

ggsave("4_manuscript/Supplementary_note/fpa_ms2/CE25/M274T390_POS.pdf", plot = p, width = 8, height = 6, dpi = 300)


# for M292T390_POS
which(ms2_data@variable_id == "M292T390_POS")

spec_1 <- ms2_data@ms2_spectra[[2325]] %>% as.data.frame()

spec_1_standard <- mpsnyder_rplc_ms2@spectra.data[["Spectra.positive"]][["RPLC_259"]][["NCE25"]]


# normalization
spec_1 <- spec_1 %>%
  mutate(rel_intensity = intensity / max(intensity, na.rm = TRUE))
spec_1_standard <- spec_1_standard %>%
  mutate(rel_intensity = intensity / max(intensity, na.rm = TRUE))

# add labels
spec_1$source <- "Sample"
spec_1_standard$source <- "Standard"
spec_1_standard$rel_intensity <- -spec_1_standard$rel_intensity

# combine
spec_all <- bind_rows(spec_1, spec_1_standard)

## calculate axis limits
mz_min <- floor(min(spec_all$mz, na.rm = TRUE))
mz_max <- ceiling(max(spec_all$mz, na.rm = TRUE))
x_pad  <- (mz_max - mz_min) * 0.02
x_lim  <- c(mz_min - x_pad, mz_max + x_pad)

p <- ggplot(spec_all, aes(x = mz, y = rel_intensity, color = source)) +  

  geom_segment(aes(xend = mz, yend = 0), linewidth = 0.53) +  # 0.35 * 1.5 ≈ 0.53
  scale_color_manual(values = c("Standard" = "red", "Sample" = "black")) +
  guides(color = "none") +
  labs(
    x = "Mass to charge ratio (m/z)",
    y = "Relative intensity"
  ) +

  scale_y_continuous(
    limits = c(-1.1, 1.1),  
    breaks = seq(-1, 1, 0.5),
    labels = function(y) abs(y)
  ) +
  coord_cartesian(xlim = x_lim, expand = FALSE) +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.ticks.length = unit(2.5, "pt"),
    plot.title = element_text(hjust = 0, face = "bold"),
    plot.subtitle = element_text(hjust = 0),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.25),
    panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60")

print(p)

ggsave("4_manuscript/Supplementary_note/fpa_ms2/CE25/M292T390_POS.pdf", plot = p, width = 8, height = 6, dpi = 300)


# for M232T328_2_POS
which(ms2_data@variable_id == "M232T328_2_POS")

spec_1 <- ms2_data@ms2_spectra[[1526]] %>% as.data.frame()

spec_1_standard <- mpsnyder_rplc_ms2@spectra.data[["Spectra.positive"]][["RPLC_399"]][["NCE25"]]


# normalization
spec_1 <- spec_1 %>%
  mutate(rel_intensity = intensity / max(intensity, na.rm = TRUE))
spec_1_standard <- spec_1_standard %>%
  mutate(rel_intensity = intensity / max(intensity, na.rm = TRUE))

# add labels
spec_1$source <- "Sample"
spec_1_standard$source <- "Standard"
spec_1_standard$rel_intensity <- -spec_1_standard$rel_intensity

# combine
spec_all <- bind_rows(spec_1, spec_1_standard)

## calculate axis limits
mz_min <- floor(min(spec_all$mz, na.rm = TRUE))
mz_max <- ceiling(max(spec_all$mz, na.rm = TRUE))
x_pad  <- (mz_max - mz_min) * 0.02
x_lim  <- c(mz_min - x_pad, mz_max + x_pad)

p <- ggplot(spec_all, aes(x = mz, y = rel_intensity, color = source)) +  
  
  geom_segment(aes(xend = mz, yend = 0), linewidth = 0.53) +  # 0.35 * 1.5 ≈ 0.53
  scale_color_manual(values = c("Standard" = "red", "Sample" = "black")) +
  guides(color = "none") +
  labs(
    x = "Mass to charge ratio (m/z)",
    y = "Relative intensity"
  ) +
  
  scale_y_continuous(
    limits = c(-1.1, 1.1),  
    breaks = seq(-1, 1, 0.5),
    labels = function(y) abs(y)
  ) +
  coord_cartesian(xlim = x_lim, expand = FALSE) +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.ticks.length = unit(2.5, "pt"),
    plot.title = element_text(hjust = 0, face = "bold"),
    plot.subtitle = element_text(hjust = 0),
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.25),
    panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60")

print(p)

ggsave("4_manuscript/Supplementary_note/fpa_ms2/CE25/M232T328_2_POS.pdf", plot = p, width = 8, height = 6, dpi = 300)

