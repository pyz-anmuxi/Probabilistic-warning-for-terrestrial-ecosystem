library(ggplot2)
library(maps)
library(readxl)
library(dplyr)
library(ggnewscale)

filename <- "D:/CODE/Exit time/Ecosystem/Data_test/Multi_data/exit_sites.xlsx"


# 读取 Excel 文件（默认读取第一个 sheet）
data <- read_excel(filename, sheet = 1)

left_larger <- filter(data, `left_basin` > `right_basin`)
right_larger <- filter(data, `left_basin` < `right_basin`)

# 获取世界地图数据 
world_map <- map_data("world")

windows(width = 10, height = 5)
p <- ggplot() +
  # 绘制世界地图背景
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "grey95", color = "grey50", linewidth = 0.7) +
  # 第一层：Left > Right 的点（绿色渐变显示left basin值）
  geom_point(data = right_larger, aes(x = long, y = lat, fill = `left_basin`),
             shape = 21, color = "black", size = 3, stroke = 0.8) +
  
  scale_fill_gradient(low = "lightgreen", high = "darkgreen",
                      name = expression(italic("T")[av]*" in high > "*italic("T")[av]*" in low"),
                      guide = guide_colorbar(order = 1)) +
  
  new_scale_fill() +
  
  # 第二层：Right ≥ Left 的点（红色渐变显示right basin值）
  geom_point(data = left_larger, aes(x = long, y = lat, fill = `right_basin`),
             shape = 21, color = "black", size = 3, stroke = 0.8) +
  scale_fill_gradient(low = "pink", high = "darkred",
                      name = expression(italic("T")[av]*" in high < "*italic("T")[av]*" in low"),
                      guide = guide_colorbar(order = 2)) +
  
  
  # 设置坐标轴格式
  scale_x_continuous(
    breaks = c(-100, 0, 100), 
    labels = c("100°W", "0°", "100°E")  # 显示 `°W`, `°E`
  ) +
  scale_y_continuous(
    breaks = c(-60, 0, 60), 
    labels = c("60°S", "0°", "60°N")  # 显示 `°S`, `°N`
  ) +
  
  # **确保地图比例正确**
  coord_fixed(ratio = 1.2, xlim = c(-180, 180), ylim = c(-60, 80)) +  
  
  theme_minimal(base_size = 14) +  # 设置基础字体大小
  theme(
    panel.border = element_rect(color = "grey50", fill = NA, linewidth = 1.2),
    axis.text = element_text(size = 12),  # 坐标轴刻度字体大小
    axis.title = element_text(size = 14),  # 坐标轴标题字体
    panel.grid.major = element_line(linewidth = 0.5, color = "gray80"), 
    legend.position = "bottom",  # 图例位置在上方
    legend.box = "horizontal",  # 水平排列多个图例
    legend.box.just = "center",  # 居中对齐
  ) +
  labs(x = "Longitude", y = "Latitude")

print(p)


####################################

# 获取世界地图数据 
world_map <- map_data("world")

windows(width = 10, height = 5)
p <- ggplot() +
  # 绘制世界地图背景
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "#f5f5f5", color = "#cccccc", size = 0.7) +
  # 第一层：Left > Right 的点（绿色渐变显示left basin值）
  geom_point(data = data, aes(x = long, y = lat, fill = `left_right`),
             shape = 21, color = "black", size = 3, stroke = 0.8) +
  
  scale_fill_gradient2(low = "#4575B4", mid = "#FFFFBF", high = "#D73027",
                       midpoint = 1,  # 设置中间值，可根据数据调整
                       expression(italic("T")[av]^{low} / italic("T")[av]^{high}),
                       guide = guide_colorbar(order = 1)) +

  # 设置坐标轴格式
  scale_x_continuous(
    breaks = c(-100, 0, 100), 
    labels = c("100°W", "0°", "100°E")  # 显示 `°W`, `°E`
  ) +
  scale_y_continuous(
    breaks = c(-60, 0, 60), 
    labels = c("60°S", "0°", "60°N")  # 显示 `°S`, `°N`
  ) +
  
  # **确保地图比例正确**
  coord_fixed(ratio = 1.2, xlim = c(-180, 180), ylim = c(-60, 80)) +  
  
  theme_minimal(base_size = 14) +  # 设置基础字体大小
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#e6f3ff", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "grey50", fill = NA, linewidth = 1.2),
    axis.text = element_text(size = 12, family = "serif"),  # 坐标轴刻度字体大小
    axis.title = element_text(size = 14, family = "serif"),  # 坐标轴标题字体
    legend.position = "bottom",  # 图例位置在上方
    legend.box = "horizontal",  # 水平排列多个图例
    legend.box.just = "center",  # 居中对齐
  ) +
  labs(x = "Longitude", y = "Latitude")

print(p)
