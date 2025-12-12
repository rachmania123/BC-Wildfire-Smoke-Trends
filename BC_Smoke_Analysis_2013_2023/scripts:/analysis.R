# ==============================================================================
# PROJECT: BC Wildfire Smoke Trends Analysis (2013-2023)
# DESCRIPTION: 
#   This script performs statistical trend analysis, satellite validation,
#   and spatial clustering analysis on wildfire smoke data.
#
# INPUTS:  data/processed_data_full.csv
# OUTPUTS: outputs/ (Figures)
# ==============================================================================

# LIBRARIES

library(tidyverse)
library(lubridate)
library(purrr)
library(sf)
library(spdep)       # For Spatial Autocorrelation
library(ggspatial)   # For Map tiles
library(viridis)     # Color palettes
library(scales)
library(patchwork)
library(Kendall)     # Mann-Kendall Trend Test
library(rstatix)     # Effect Size
library(pROC)        # ROC Curve analysis
library(stringr)

# Set plotting theme
theme_set(theme_minimal(base_size = 14) +
            theme(
              panel.grid.minor = element_blank(),
              plot.title = element_text(face = "bold"),
              strip.text = element_text(face = "bold", size = 12)
            ))



# LOAD DATAFRAME

df <- read_csv("data/processed_data_full.csv", show_col_types = FALSE)

# Clean and Classify for Analysis
df_main <- df %>%
  filter(!is.na(smoke_PM25), !is.na(PM25)) %>%
  mutate(
    date = make_date(Year, Month, Day),
    
    # 1. Classify Area Type
    Area_Type = case_when(
      str_detect(City, "Vancouver|Burnaby|Richmond|Surrey") | City %in% c("Victoria", "Kelowna") ~ "Urban",
      str_detect(City, "Metro Van") | City %in% c("Abbotsford", "Kamloops", "Nanaimo", "Prince George", 
                                                  "Chilliwack", "Vernon", "Courtenay", "Campbell River", 
                                                  "Penticton", "Mission", "Maple Ridge", "New Westminster") ~ "Suburban",
      TRUE ~ "Rural"
    ),
    Area_Type = factor(Area_Type, levels = c("Urban", "Suburban", "Rural")),
    
    # 2. Define Threshold for Smoke Days (> 25 ug/m3)
    Is_Unhealthy = ifelse(smoke_PM25 > 25, 1, 0)
  )



# RQ1: Variation Across Years 
kw_year <- kruskal.test(smoke_PM25 ~ Year, data = df_main)
print("Kruskal-Wallis Test Results (Years)")
print(kw_year)

# RQ2: Trends (Mann-Kendall & Linear Slope) 

# A. Intensity Trend (Mean PM2.5)
annual_avg <- df_main %>% group_by(Year) %>% summarise(Mean_Smoke = mean(smoke_PM25, na.rm=TRUE))
mk_intensity <- MannKendall(annual_avg$Mean_Smoke)
print("Mann-Kendall (Intensity)")
summary(mk_intensity)

# B. Frequency Trend (Smoke Days)
trend_data <- df_main %>%
  filter(HMS_any == 1) %>% 
  group_by(Year) %>%
  summarise(Smoke_Days = n_distinct(date), .groups = "drop") %>%
  arrange(Year)

mk_proxy <- cor.test(trend_data$Year, trend_data$Smoke_Days, method = "kendall")
lm_model <- lm(Smoke_Days ~ Year, data = trend_data)
slope_val <- coef(lm_model)[["Year"]]

cat("Kendall's Tau (Freq):", mk_proxy$estimate, "\n")
cat("P-value (Freq):", mk_proxy$p.value, "\n")
cat("Estimated Slope (Days/Year):", slope_val, "\n")

#  RQ3: Satellite Validation Interaction Model 
print(" Interaction Test (Satellite Accuracy by Area) ")
interaction_model <- lm(smoke_PM25 ~ MAIAC_AOD_550 * Area_Type, data = df_main)
print(anova(interaction_model))

print(paste("Overall Correlation:", cor(df_main$MAIAC_AOD_550, df_main$smoke_PM25, use="complete.obs")))

#  RQ4: Urban vs Rural Comparison 
kw_area <- kruskal.test(smoke_PM25 ~ Area_Type, data = df_main)
print(" Kruskal-Wallis Test Results (Urban vs Rural) ")
print(kw_area)

if(kw_area$p.value < 0.05) {
  print(" Post-Hoc Pairwise Wilcoxon ")
  print(pairwise.wilcox.test(df_main$smoke_PM25, df_main$Area_Type, p.adjust.method = "bonferroni"))
}

# Effect Size
effect_size_area <- kruskal_effsize(smoke_PM25 ~ Area_Type, data = df_main)
print(" Effect Size (Urban vs Rural) ")
print(effect_size_area)

# Spatial Autocorrelation (Moran's I) 

# Prepare Spatial Data
station_spatial <- df_main %>%
  mutate(geo_content = str_extract(.geo, "(?<=\\[).+?(?=\\])")) %>%
  separate(geo_content, into = c("lon_str", "lat_str"), sep = ",", remove = FALSE) %>%
  mutate(lon = as.numeric(lon_str), lat = as.numeric(lat_str)) %>%
  filter(!is.na(lon), !is.na(lat)) %>%
  group_by(Station, lat, lon) %>%
  summarise(Mean_Smoke = mean(smoke_PM25, na.rm = TRUE), .groups = "drop") %>%
  distinct(lat, lon, .keep_all = TRUE) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Neighbors & Weights
knn <- knearneigh(st_coordinates(station_spatial), k = 5) 
listw <- nb2listw(knn2nb(knn), style = "W")
moran_test <- moran.test(station_spatial$Mean_Smoke, listw)
print(moran_test)



# 6. VISUALIZATIONS


#  FIGURE 1: Trend Analysis 
df_trend_panel <- df_main %>%
  filter(!is.na(background_PM25)) %>%
  group_by(Year) %>%
  summarise(
    Conc_ugm3 = mean(smoke_PM25, na.rm=TRUE),
    Contrib_Pct = (sum(smoke_PM25)/sum(PM25))*100,
    Exceed_Pct = (sum(PM25 > 25 & background_PM25 <= 25)/n())*100
  ) %>%
  pivot_longer(cols = -Year, names_to = "Metric_Type", values_to = "Value") %>%
  mutate(
    Metric_Label = case_when(
      Metric_Type == "Conc_ugm3" ~ "Average Smoke PM2.5 Concentration (µg/m³)",
      Metric_Type == "Contrib_Pct" ~ "% Smoke Contribution to Total PM2.5",
      Metric_Type == "Exceed_Pct" ~ "% Days Smoke Pushes Total PM2.5 > 25 µg/m³"
    ),
    Metric_Label = factor(Metric_Label, levels = c("Average Smoke PM2.5 Concentration (µg/m³)",
                                                   "% Smoke Contribution to Total PM2.5",
                                                   "% Days Smoke Pushes Total PM2.5 > 25 µg/m³")),
    Label_Text = ifelse(Metric_Type == "Conc_ugm3", sprintf("%.1f", Value), sprintf("%.1f%%", Value))
  )

p_3panel <- ggplot(df_trend_panel, aes(x = Year, y = Value)) +
  geom_line(color = "#D71920", linewidth = 1.2) +
  geom_point(color = "#D71920", size = 4) +
  geom_text(aes(label = Label_Text), vjust = -1.2, fontface = "bold", size = 3.5) +
  facet_wrap(~Metric_Label, ncol = 1, scales = "free_y", strip.position = "top") +
  scale_x_continuous(breaks = 2013:2023) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.3))) +
  labs(title = "Smoke PM2.5 Impact (2013-2023)") +
  theme(axis.title = element_blank())

print(p_3panel)



#  FIGURE 2: Time Series Error Bars 
station_annual <- df_main %>%
  group_by(Station, Area_Type, Year) %>%
  summarise(
    Total_Smoke_Mass = sum(smoke_PM25, na.rm = TRUE),
    Total_PM_Mass = sum(PM25, na.rm = TRUE),
    Pct_Contribution = (Total_Smoke_Mass / Total_PM_Mass) * 100, .groups = "drop"
  )

summary_ts <- station_annual %>%
  group_by(Year, Area_Type) %>%
  summarise(
    Mean_Pct = mean(Pct_Contribution, na.rm = TRUE),
    SD_Pct = sd(Pct_Contribution, na.rm = TRUE),
    Lower = pmax(Mean_Pct - SD_Pct, 0), 
    Upper = pmin(Mean_Pct + SD_Pct, 100), .groups = "drop"
  )

p_errorbar <- ggplot(summary_ts, aes(x = Year, y = Mean_Pct, color = Area_Type, group = Area_Type)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, alpha = 0.5, position = position_dodge(0.3)) +
  geom_line(size = 1, position = position_dodge(0.3)) +
  geom_point(size = 3, position = position_dodge(0.3)) +
  scale_color_manual(values = c("Urban" = "#440154", "Suburban" = "#21918c", "Rural" = "#fde725")) +
  scale_x_continuous(breaks = 2013:2023) +
  labs(title = "Annual % Smoke Contribution to Total PM2.5", y = "Smoke Contribution (%)", color = "Area Type") +
  theme(legend.position = "bottom")

print(p_errorbar)



#  FIGURE 3: Hexagon Density (Satellite Correlation) 
valid_pairs <- df_main %>% filter(!is.na(smoke_PM25), !is.na(MAIAC_AOD_550))

p_hex <- ggplot(valid_pairs, aes(x = MAIAC_AOD_550, y = smoke_PM25)) +
  geom_hex(bins = 50) +
  geom_smooth(method = "lm", color = "white", linetype = "dashed", size = 1) +
  scale_fill_viridis_c(option = "magma", name = "Count") +
  coord_cartesian(ylim = c(0, 150), xlim = c(0, 1.5)) +
  labs(title = "RQ3: Satellite AOD vs Smoke PM2.5 (Density)", x = "MAIAC AOD", y = "Smoke PM2.5 (µg/m³)")

print(p_hex)



#  FIGURE 4: Violin Plot (Distribution) 
stats_labels <- df_main %>%
  group_by(Area_Type) %>%
  summarise(
    Median = median(smoke_PM25, na.rm = TRUE),
    Q1 = quantile(smoke_PM25, 0.25, na.rm = TRUE),
    Q3 = quantile(smoke_PM25, 0.75, na.rm = TRUE),
    Label = paste0("Med: ", round(Median, 1), "\n(", round(Q1, 1), "-", round(Q3, 1), ")"),
    y_pos = 18
  )

p_violin <- ggplot(df_main, aes(x = Area_Type, y = smoke_PM25, fill = Area_Type)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.8, outlier.shape = NA) +
  geom_text(data = stats_labels, aes(y = y_pos, label = Label), color = "black", size = 3.5, fontface = "bold", vjust = 0) +
  scale_fill_viridis_d(option = "mako", begin = 0.3, end = 0.8) +
  coord_cartesian(ylim = c(0, 20)) +
  labs(title = "RQ4: Smoke Exposure Distribution", x = NULL, y = "Smoke PM2.5 (µg/m³)") +
  theme(legend.position = "none")

print(p_violin)



#  FIGURE 5: Spatial Risk Map 
station_map_data <- df_main %>%
  mutate(geo_content = str_extract(.geo, "(?<=\\[).+?(?=\\])")) %>%
  separate(geo_content, into = c("lon_str", "lat_str"), sep = ",", remove = FALSE) %>%
  mutate(lon = as.numeric(lon_str), lat = as.numeric(lat_str)) %>%
  filter(!is.na(lon), !is.na(lat)) %>%
  group_by(Station, Area_Type, lat, lon) %>%
  summarise(Pct_Unhealthy = mean(Is_Unhealthy) * 100, .groups = "drop") %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(3857)

bbox <- st_bbox(station_map_data)

p_map <- ggplot() +
  annotation_map_tile(type = "cartolight", zoomin = 0, alpha = 0.8) +
  geom_sf(data = station_map_data, aes(fill = Pct_Unhealthy, shape = Area_Type), size = 5, color = "black", stroke = 0.5) +
  scale_fill_gradientn(colors = c("white", "#fed976", "#fd8d3c", "#e31a1c", "#800026"), name = "% Days > 25µg/m³") +
  scale_shape_manual(values = c("Urban" = 22, "Suburban" = 24, "Rural" = 21), name="Area Type") +
  coord_sf(crs = 3857, xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tl", style = north_arrow_minimal) +
  labs(title = "Spatial Distribution of Smoke Risk in BC", x = NULL, y = NULL) +
  theme(axis.text = element_blank())

print(p_map)



#  FIGURE 6: ROC Curve 
roc_data <- valid_pairs %>% mutate(High_Smoke = ifelse(smoke_PM25 > 25, 1, 0))
roc_obj <- roc(roc_data$High_Smoke, roc_data$MAIAC_AOD_550, quiet = TRUE)

p_roc <- ggroc(roc_obj) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "red") +
  labs(title = paste0("RQ5: ROC Curve (AUC = ", round(auc(roc_obj), 3), ")"))

print(p_roc)



#  FIGURE 7: RQ2 Trend Plot (Smoke Days) 
p_trend_days <- ggplot(trend_data, aes(x = Year, y = Smoke_Days)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", color = "red", linetype = "dashed") +
  labs(title = paste0("RQ2: Trend in Smoke Days (tau=", round(mk_proxy$estimate, 2), ")"),
       subtitle = paste("Slope ~", round(slope_val, 2), "days/year"))

print(p_trend_days)

#  FIGURE 8: RQ1 Variation Zoomed 
p_var_zoom <- ggplot(df_main, aes(x = factor(Year), y = smoke_PM25)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue", alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
  coord_cartesian(ylim = c(0, 10)) +
  labs(title = "RQ1: Annual Variation in Smoke PM2.5 (Zoomed)", x = "Year", y = "Smoke PM2.5 (µg/m³)")

print(p_var_zoom)
