library(ncdf4)
library(ggplot2)

# ── 1. Read the NetCDF ─────────────────────────────────────────────────────────
nc  <- nc_open("landcoverModisGlobal_2020_1D.nc")
lc_s  <- ncvar_get(nc, "landcover")   # [lon × lat]
lon_s <- ncvar_get(nc, "lon")
lat_s <- ncvar_get(nc, "lat")
nc_close(nc)

# ── 2. MODIS MCD12Q1 Type 1 class labels (0–17) ────────────────────────────────
class_labels <- c(
  "0"  = "Water",
  "1"  = "Evergreen Needleleaf Forest",
  "2"  = "Evergreen Broadleaf Forest",
  "3"  = "Deciduous Needleleaf Forest",
  "4"  = "Deciduous Broadleaf Forest",
  "5"  = "Mixed Forest",
  "6"  = "Closed Shrublands",
  "7"  = "Open Shrublands",
  "8"  = "Woody Savannas",
  "9"  = "Savannas",
  "10" = "Grasslands",
  "11" = "Permanent Wetlands",
  "12" = "Croplands",
  "13" = "Urban & Built-up",
  "14" = "Cropland / Natural Mosaic",
  "15" = "Snow & Ice",
  "16" = "Barren",
  "17" = "Unclassified"
)

# Palette: one distinct colour per class
class_colors <- c(
  "Water"                      = "#FFFFFF",
  "Evergreen Needleleaf Forest"= "#1A5C38",
  "Evergreen Broadleaf Forest" = "#2D8B4E",
  "Deciduous Needleleaf Forest"= "#57A058",
  "Deciduous Broadleaf Forest" = "#8CC56A",
  "Mixed Forest"               = "#B2DF6A",
  "Closed Shrublands"          = "#C29A3F",
  "Open Shrublands"            = "#DDB96E",
  "Woody Savannas"             = "#A8D29C",
  "Savannas"                   = "#D0E8A0",
  "Grasslands"                 = "#F5F5A0",
  "Permanent Wetlands"         = "#7FC4C4",
  "Croplands"                  = "#F2C96A",
  "Urban & Built-up"           = "#E84040",
  "Cropland / Natural Mosaic"  = "#F0B030",
  "Snow & Ice"                 = "#DDDDDD",
  "Barren"                     = "#C0A080",
  "Unclassified"               = "#FFFFFF"
)

# ── 3. Build data frame ────────────────────────────────────────────────────────
# lc matrix is [lon(3600) × lat(1494)]; we need [lat × lon] for image plotting

# Melt to long format
df <- data.frame(
  lon   = rep(lon_s, times = length(lat_s)),          # lon varies fastest
  lat   = rep(lat_s, each  = length(lon_s)),           # lat is outer
  class = as.integer(as.vector(lc_s))                  # row-major: lon inner
)
df <- df[!is.na(df$class), ]
df$label <- factor(class_labels[as.character(df$class)],
                   levels = names(class_colors))

# ── 4. World outline via maps ──────────────────────────────────────────────────
library(maps)
world <- map_data("world")

# ── 5. Plot ────────────────────────────────────────────────────────────────────
pal <- setNames(class_colors, unname(class_labels))
windows(width = 10, height = 5)
ggplot() +
  geom_raster(data = df, aes(x = lon, y = lat, fill = label)) +
  geom_path(data = world,
            aes(x = long, y = lat, group = group),
            colour = "grey20", linewidth = 0.2) +
  scale_fill_manual(
    values   = class_colors,
    na.value = "transparent",
    name     = "Land Cover",
    guide    = guide_legend(ncol = 1, keywidth = 0.8, keyheight = 0.65,
                            override.aes = list(size = 3))
  ) +
  coord_fixed(1.3, xlim = c(-180, 180), ylim = c(-60, 80)) +
  labs(
    title    = "MODIS Global Land Cover 2020",
    subtitle = "MCD12Q1 · IGBP Classification",
    x = "Longitude", y = "Latitude",
    caption  = "Source: NASA MODIS / Esri ArcGIS"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, colour = "grey40"),
    legend.position = "right",
    legend.text   = element_text(size = 7),
    legend.title  = element_text(size = 9, face = "bold"),
    panel.grid    = element_blank(),
    axis.text     = element_text(size = 8)
  )


################
library(ggplot2)
library(maps)
library(cowplot)
library(ggnewscale)
filename <- "D:/CODE/Exit time/Ecosystem/Data_test/Multi_data/exit_sites.xlsx"


# 读取 Excel 文件（默认读取第一个 sheet）
data <- read_excel(filename, sheet = 1)

# 获取世界地图数据 
world_map <- map_data("world")

windows(width = 10, height = 5)

p_main <- ggplot() +
  # 第一层：植被类型
  geom_raster(data = df, aes(x = lon, y = lat, fill = label)) +
  scale_fill_manual(values = class_colors, name = "Land Cover",
                    guide = guide_legend(
                      ncol = 1, 
                      keywidth = 0.8, 
                      keyheight = 0.65,
                      title.position = "top", 
                      label.position = "right")) +
  
  new_scale_fill() +
  # 第二层：世界地图边界
  geom_path(data = world_map,
            aes(x = long, y = lat, group = group),
            colour = "grey30", linewidth = 0.3) +
  # 第三层：
  geom_point(data = data, aes(x = long, y = lat, fill = `left_right`),
             shape = 24,color = "black",size = 1.5, stroke = 0.8) +
  
  
  scale_fill_gradient2(low = "#4575B4", mid = "#FFFFBF", high = "#D73027",
                       midpoint = 1,
                       expression(italic("T")[av]^{low} / italic("T")[av]^{high}),
                       guide = guide_colorbar(
                         title.position = "top",
                         barwidth = 20,
                         barheight = 0.8,
                         label.position = "bottom"
                       )) +
  # 设置坐标轴格式
  scale_x_continuous(
    breaks = c(-100, 0, 100), 
    labels = c("100°W", "0°", "100°E")  # 显示 `°W`, `°E`
  ) +
  scale_y_continuous(
    breaks = c(-60, 0, 60), 
    labels = c("60°S", "0°", "60°N")  # 显示 `°S`, `°N`
  ) +
  coord_fixed(ratio = 1.2, xlim = c(-180, 180), ylim = c(-60, 80)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",  # 移除主图的所有图例
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#FFFFFF", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "grey50", fill = NA, linewidth = 1.2),
    axis.text = element_text(size = 12, family = "serif"),  # serif 通常映射到 Times New Roman
    axis.title = element_text(size = 14, family = "serif")
  ) +
  labs(x = "Longitude", y = "Latitude")

# 2. 提取颜色渐变图例（底部居中）
legend_color <- get_legend(
  ggplot() +
    geom_point(data = data, aes(x = long, y = lat, color = `left_right`)) +
    scale_color_gradient2(
      low = "#4575B4", mid = "#FFFFBF", high = "#D73027",
      midpoint = 1,
      name = expression(italic("T")[av]^{low} / italic("T")[av]^{high}),
      guide = guide_colorbar(
        title.position = "left",
        barwidth = unit(8, "cm"),
        barheight = unit(0.8, "cm"),
        label.position = "bottom"
      )
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      plot.margin = margin(0, 0, 0, 0),
      text = element_text(family = "serif")
    )
)

# 3. 提取植被类型图例（右侧）
legend_fill <- get_legend(
  ggplot() +
    geom_raster(data = df, aes(x = lon, y = lat, fill = label)) +
    scale_fill_manual(
      values = class_colors,
      name = "Land Cover",
      guide = guide_legend(
        ncol = 1,
        keywidth = 0.8,
        keyheight = 0.65,
        title.position = "top",
        label.position = "right"
      )
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.justification = "center",
      plot.margin = margin(0, 0, 0, 10),
      text = element_text(family = "serif")
    )
)

# 4. 组合图形
# 先将主图和右侧图例组合
p_with_right <- plot_grid(
  p_main,
  legend_fill,
  ncol = 2,
  rel_widths = c(1, 0.25)  # 右侧图例占25%宽度
)

# 再在底部添加颜色图例（居中）
p_final <- plot_grid(
  p_with_right,
  legend_color,
  ncol = 1,
  rel_heights = c(1, 0.12)  # 底部图例占12%高度
)

print(p_final)

