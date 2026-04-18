library(ggplot2)
library(maps)
library(dplyr)

data <- read.table("./0_Gimms_kNDVI_validated_sites.txt",head=FALSE)
#data <- read.table("./pku_Gimms_kNDVI_validated_sites.txt",head=FALSE)
longitude <- data[2,]
latitude <- data[1,]
values <- data[3,]

df <- data.frame(
  lon = as.numeric(longitude),
  lat = as.numeric(latitude),
  value = as.numeric(values)
)
df_red <- df[df$value == 1, ]    # 红色：值为1
df_blue <- df[df$value == -1, ]  # 蓝色：值为-1

locations <- data.frame(
  lat = c(49.25, 16.25), 
  lon = c(138.25, -95.75), 
  name = c("health", "un")
)


world_map <- map_data("world")
windows(width = 10, height = 5)
ggplot() +
  # 绘制世界地图底图
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "#f5f5f5", color = "#cccccc", size = 0.7) +
  
  # 绘制值为1的点（红色）
  geom_point(data = df, aes(x = lon, y = lat, fill = factor(value)),
             shape = 21, color = "black", size = 3, stroke = 0.8) +
  scale_fill_manual(values = c("1" = "orange", "-1" = "orange"))+

  geom_point(data = locations[1,], aes(x = lon, y = lat, fill = name),
             shape = 24, fill="#CE0202", color = "black", size = 4, stroke = 0.8) + 
  
  geom_point(data = locations[2,], aes(x = lon, y = lat, fill = name),
             shape = 24, fill="#65b168", color = "black", size = 4, stroke = 0.8) + 
  
  # 设置坐标轴格式
  scale_x_continuous(
    breaks = c(-100, 0, 100), 
    labels = c("100°W", "0°", "100°E")  # 显示 `°W`, `°E`
  ) +
  scale_y_continuous(
    breaks = c(-60, 0, 60), 
    labels = c("60°S", "0°", "60°N")  # 显示 `°S`, `°N`
  ) +
  # 设置坐标范围（根据你的数据调整）
  coord_fixed(ratio = 1.2, xlim = c(-180, 180), ylim = c(-60, 80)) +
  # 主题设置
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#e6f3ff", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),  # 坐标轴刻度字体大小
    axis.title = element_text(size = 14),  # 坐标轴标题字体
    panel.border = element_rect(color = "grey50", fill = NA, linewidth = 1.2)
  ) +
  labs(x = "Longitude", y = "Latitude")





