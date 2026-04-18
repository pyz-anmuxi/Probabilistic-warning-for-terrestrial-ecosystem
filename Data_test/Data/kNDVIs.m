clc;clear;close all;

file = 'D:\CODE\Exit time\Ecosystem\Data_test\Data\';
FilePath = '0_Gimms_kNDVI_validated_sites.txt';

var = load(strcat(file,FilePath));
%%
close all;

tab = var(4:end,1:1);

tim = 1: 1: size(tab);
xi = 1: 0.1: size(tab);
tab = interp1(tim, tab, xi, 'PCHIP'); % 线性插值

namex = 'Vegetation state, V';
timeunit = '(Month)';

plot(xi, tab, 'Color',[170 170 170]/255,'LineWidth',2)

ylim([min(tab)-0.05 max(tab)+0.05]);
set(gca, 'XTick',[3*12 13*12 23*12 33*12], 'XTickLabel', {'1985','1995','2005','2015'});
xlim([tim(1) tim(end)]);
set(gcf, 'position', [200 200 900 450]);
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
xlabel(['Time ' timeunit],'fontsize',18,'fontname','Times')	%设置y轴字体
ylabel(namex,'fontsize',18,'fontname','Times')	%设置x轴字体
box on


