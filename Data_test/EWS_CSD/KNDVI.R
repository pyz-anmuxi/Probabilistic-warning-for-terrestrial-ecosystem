library(ncdf4)
library(EWSmethods)

ncfile <- nc_open("D:/CODE/Exit time/Ecosystem/Data_test/Data/kNDVI_Example.nc")
data <- ncvar_get(ncfile, "var")


kNDVI_month <- data[5:400,4:4]

#kNDVI_month = (tab - mean(tab)) / sd(tab)
plot(kNDVI_month, type = "l", lwd = 2) 

time <- seq(1,length(kNDVI_month))
NDVI_phase1 <- data.frame(time = time, kNDVI = kNDVI_month)
rolling_ews_eg <- uniEWS(data = NDVI_phase1,
                          metrics = c("ar1","SD","cv","skew","kurt","rr"),
                          method = "rolling",
                          winsize = 5)
plot(rolling_ews_eg,  y_lab = "Density")

EWS <- rolling_ews_eg[["EWS"]]

ews <- EWS[["raw"]]

# 保存为 .txt 文件（制表符分隔）
write.table(ews, file = "ndvi_2.txt", sep = "\t", row.names = FALSE)





library(zoo)
library(changepoint)

min_period <- 72
cpt_result <- cpt.var(kNDVI_month,
                      method = "PELT",
                      penalty = "MBIC",
                      Q = 2,
                      minseglen = min_period)
plot(cpt_result)
# 排除首尾72个月内的变点
detected_cpts <- cpts(cpt_result)
# 排除首尾72个月内的变点
valid_cpts <- detected_cpts[detected_cpts > min_period & 
                              detected_cpts < (length(kNDVI_month) - min_period)]

tipping_point <- cpt_result@cpts  # 变点结果

phase1 <- kNDVI_month[1:tipping_point[1]]
plot(phase1, type = "l", lwd = 2)

time <- seq(1,length(phase1))
NDVI_phase1 <- data.frame(time = time, kNDVI = phase1)
rolling_ews_eg1 <- uniEWS(data = NDVI_phase1,
                         metrics = c("ar1","SD","cv","skew","kurt","rr"),
                         method = "rolling",
                         winsize = 5)
x(11)  
plot(rolling_ews_eg1,  y_lab = "Density")


phase2 <- kNDVI_month[tipping_point[1]:tipping_point[2]]
plot(phase2, type = "l", lwd = 2)

time <- seq(1,length(phase2))
NDVI_phase2 <- data.frame(time = time, kNDVI = phase2)
rolling_ews_eg2 <- uniEWS(data = NDVI_phase2,
                          metrics = c("ar1","SD","skew"),
                          method = "rolling",
                          winsize = 5)
x(11)  
plot(rolling_ews_eg2,  y_lab = "Density")
