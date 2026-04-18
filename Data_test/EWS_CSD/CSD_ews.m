clc; clear;close all

E1 = readtable('D:\CODE\Exit time\Ecosystem\EWS_CSD\ndvi_1.txt', 'Delimiter', '\t');  % 制表符分隔
E2 = readtable('D:\CODE\Exit time\Ecosystem\EWS_CSD\ndvi_2.txt', 'Delimiter', '\t');  % 制表符分隔
%%  ndvi_1

ndvi = E1.count_used;
ndvi = (ndvi - mean(ndvi)) / std(ndvi);
time = E1.time;
plot(time, ndvi, 'Color',[130/255 177/255 255/255],'LineWidth',1.2)
line([1, 150], [0.047,0.047], 'LineWidth', 2, 'LineStyle', '-.', 'color', [0.6350, 0.0780, 0.1840]);
line([150, 200], [1.8, 1.8], 'LineWidth', 2, 'LineStyle', '-.', 'color',  [55/256 131/256 59/256]);
line([200, 396], [0.047, 0.047],  'LineWidth', 2, 'LineStyle', '-.', 'color', [0.6350, 0.0780, 0.1840]);

ylim([min(ndvi)-0.05, max(ndvi)+0.05]);
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
ylabel('Vegetation state, V','fontsize',18,'fontname','Times')	%设置x轴字体
set(gca, 'XTick',[100 200 300 400]);  % 隐藏x轴刻度标签
box on
set(gcf, 'position', [500 100 900 250]);

%%
close all;

subplot(3,2,1)
plot(time, E1.ar1, 'Color',[165 28 54]/255,'LineWidth',1.2)
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
ylabel('AR1','fontsize',18,'fontname','Times')	%设置x轴字体
set(gca, 'XTick',[100 200 300 400], 'XTickLabel', {});  % 隐藏x轴刻度标签
box on

subplot(3,2,2)
plot(time, E1.rr, 'Color',[122 187 219]/255,'LineWidth',1.2)
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
ylabel('Return rate','fontsize',18,'fontname','Times')	%设置x轴字体
set(gca, 'XTick',[100 200 300 400], 'XTickLabel', {});  % 隐藏x轴刻度标签
box on

subplot(3,2,3)
plot(time, E1.SD, 'Color',[132 186 66]/255,'LineWidth',1.2)
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
ylabel('SD','fontsize',18,'fontname','Times')	%设置x轴字体
set(gca, 'XTickLabel', {});  % 隐藏x轴刻度标签
box on


subplot(3,2,4)
plot(time, E1.cv, 'Color',[104 36 135]/255,'LineWidth',1.2)
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
ylabel('CV','fontsize',18,'fontname','Times')	%设置x轴字体
set(gca, 'XTick',[100 200 300 400], 'XTickLabel', {});  % 隐藏x轴刻度标签
box on

%
subplot(3,2,5)
plot(time, E1.skew, 'Color',[219 180 40]/256,'LineWidth',1.2)
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
xlabel('Time','fontsize',18,'fontname','Times')	%设置x轴字体
ylabel('Skewness','fontsize',18,'fontname','Times')	%设置x轴字体
set(gca, 'XTick',[100 200 300 400]);  % 隐藏x轴刻度标签
box on

subplot(3,2,6)
plot(E1.time, E1.kurt, 'Color',[212 86 46]/255,'LineWidth',1.2)
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
xlabel('Time','fontsize',18,'fontname','Times')	%设置x轴字体
ylabel('Kurtosis','fontsize',18,'fontname','Times')	%设置x轴字体
set(gca, 'XTick',[100 200 300 400]);  % 隐藏x轴刻度标签
box on
set(gcf, 'position', [500 100 900 600]);

%%  ndvi_2
ndvi2 = E2.count_used;
ndvi2 = (ndvi2 - mean(ndvi2)) / std(ndvi2);
time = E2.time;
plot(time, ndvi2, 'Color',[130/255 177/255 255/255],'LineWidth',1.2)
line([1, 197], [0.5, 0.5], 'LineWidth', 2, 'LineStyle', '-.', 'color', [55/256 131/256 59/256]);
line([196, 298], [-1.2, -1.2], 'LineWidth', 2, 'LineStyle', '-.', 'color',  [0.6350, 0.0780, 0.1840]);
line([298, 396], [0.5, 0.5],  'LineWidth', 2, 'LineStyle', '-.', 'color', [55/256 131/256 59/256]);
ylim([min(ndvi2)-0.05, max(ndvi2)+0.05]);
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
ylabel('Vegetation state, V','fontsize',18,'fontname','Times')	%设置x轴字体
set(gca, 'XTick',[100 200 300 400]);  % 隐藏x轴刻度标签
box on
set(gcf, 'position', [500 100 900 250]);
%%
figure(2)
subplot(3,2,1)
plot(time, E2.ar1, 'Color',[165 28 54]/255,'LineWidth',1.2)
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
ylabel('AR1','fontsize',18,'fontname','Times')	%设置x轴字体
set(gca, 'XTick',[100 200 300 400], 'XTickLabel', {});  % 隐藏x轴刻度标签
box on

subplot(3,2,2)
plot(time, E2.rr, 'Color',[122 187 219]/255,'LineWidth',1.2)
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
ylabel('Return rate','fontsize',18,'fontname','Times')	%设置x轴字体
set(gca, 'XTick',[100 200 300 400], 'XTickLabel', {});  % 隐藏x轴刻度标签
box on

subplot(3,2,3)
plot(time, E2.SD, 'Color',[132 186 66]/255,'LineWidth',1.2)
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
ylabel('SD','fontsize',18,'fontname','Times')	%设置x轴字体
set(gca, 'XTick',[100 200 300 400], 'XTickLabel', {});  % 隐藏x轴刻度标签
box on

subplot(3,2,4)
plot(time, E2.cv, 'Color',[104 36 135]/255,'LineWidth',1.2)
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
ylabel('CV','fontsize',18,'fontname','Times')	%设置x轴字体
set(gca, 'XTick',[100 200 300 400], 'XTickLabel', {});  % 隐藏x轴刻度标签
box on

%
subplot(3,2,5)
plot(time, E2.skew, 'Color',[219 180 40]/256,'LineWidth',1.2)
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
xlabel('The number of data','fontsize',18,'fontname','Times')	%设置x轴字体
ylabel('Skewness','fontsize',18,'fontname','Times')	%设置x轴字体
set(gca, 'XTick',[100 200 300 400]);  % 隐藏x轴刻度标签
box on

subplot(3,2,6)
plot(time, E2.kurt, 'Color',[212 86 46]/255,'LineWidth',1.2)
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
xlabel('The number of data','fontsize',18,'fontname','Times')	%设置x轴字体
ylabel('Kurtosis','fontsize',18,'fontname','Times')	%设置x轴字体
set(gca, 'XTick',[100 200 300 400]);  % 隐藏x轴刻度标签
box on
set(gcf, 'position', [500 100 900 600]);
