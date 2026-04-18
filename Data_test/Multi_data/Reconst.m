
% Implements Langevin Reconstruction based on k-NDVI data
% 2025.4.15

clc;clear;close all;

close all;
file = 'D:\CODE\Exit time\Ecosystem\Data_test\Data\';
FilePath = '0_Gimms_kNDVI_validated_sites.txt';
% FilePath = 'pku_Gimms_kNDVI_validated_sites.txt';

var = load(strcat(file,FilePath));

% 查看不同类型数据：
type_1 = var(:,var(3:3,:)==1);
type_2 = var(:,var(3:3,:)==-1);

%%
close all;
num = 35;
tab = var(4:end,num:num);
fprintf('lat:%g\n  long : %g \n  type:%g \n', var(1:1,num:num), var(2:2,num:num), var(3:3,num:num));
% Data preprocessing: Z-score normalization
tab = (tab - mean(tab)) / std(tab);

tim = 1: 1: size(tab);
xi = 1: 0.5: size(tab);
tab = interp1(tim, tab, xi, 'PCHIP');  % Linear interpolation

plot(tab)
set(gcf, 'position', [200 200 900 450]);


%% Step2：重构Langevin方程，得到D1和D2数据
close all;
L = min(tab); 
R = max(tab); 

bins = 8;  % 决定M_1数据量，即拟合D_1曲线光滑度，如何调？ 
% 过大：回归设计矩阵秩亏，无法达到机器精度
% 过小：曲线不光滑
Tau = 1:3;  % fitlm 拟合M(i)和Tau,影响误差
% a list of time lags considered (integer values)
% (we chose 1:5. For high resolution data we can consider bigger number of time lags).
% x(t+tau) - x(t)
dt = 0.5; % time step of the data series 影响纵轴的数量级
method = 'Nadaraya-Watson'; 

res = LangevinReconst(tab, L, R, bins, Tau, dt, method, 'kNDVI');

% save('ni','res', 'tim', 'tab', '-v7.3')

% 简单查看重构效果

mod = langevin_eq(res, 'weightedspline');

mod.namex = 'kNDVI';
mod.timeunit = '(year)';

figure(2)
mod.plot('D1', 'error', true)
hold on
mod.plot('equilibria');
plot(res.C,res.D1)

% figure(3)   
% mod.plot('D2', 'error', true)

fprintf('low:%g\n  high : %g \n ', mod.equilibria(1).x, mod.equilibria(3).x);
%% 时间序列图
close all;
namex = 'Vegetation state, V';
timeunit = '(month)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 时间序列图 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(xi, tab, 'Color',[130/255 177/255 255/255],'LineWidth',1.2)

% % hold on
yline(mod.equilibria(1).x, '-',  'LineWidth', 1.8, 'color', [255/256 183/256 77/256]);
yline(mod.equilibria(2).x, '--',  'LineWidth', 1.8, 'color', [55/256 131/256 59/256]);
yline(mod.equilibria(3).x, '-', 'LineWidth', 1.8, 'color', [255/256 183/256 77/256]);

ylim([min(tab)-0.02 max(tab)+0.02]);
% set(gca, 'Ytick', 0:2500:5000)
set(gca, 'XTick',[3*12 13*12 23*12 33*12], 'XTickLabel', {'1985','1995','2005','2015'});
xlim([tim(1) tim(end)]);

% set(findobj(gcf, '-property', 'fontsize'), 'fontsize', 10);
set(gcf, 'position', [200 500 600 240]);
% text(gca, 0.9, 0.95, 'A', 'Units', 'normalized', 'fontsize', 10, 'tag','fignotext');
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
ylabel(namex,'fontsize',18,'fontname','Times')	%设置x轴字体
xlabel(['Time ' timeunit],'fontsize',18,'fontname','Times')	%设置y轴字体
box on

%% 读取D1 和 D2 数据
% nreplicates = 0; % set to zero for run without monte carlo error propagation
% if nreplicates > 0
%     results = load('Infectives_boot1');
% else
%     results = load('.\20230708\Measles.mat');
% end

% tim = results.tim;
% data = results.data;
% res = results.results;
close all;
mod = langevin_eq(res, 'weightedspline');

mod.namex = 'Vegetation state, V';
mod.timeunit = '(month)';
% 
% fprintf('Small scale is %f\n', mod.equilibria(1).x)
% fprintf('scale threshold is %f\n', mod.equilibria(2).x)
% fprintf('big scale is %f\n', mod.equilibria(3).x * (max(tab_1) - min(tab_1))+ min(tab_1))
% *************** ***************Drift: D_1 *************** ***************
subplot(3,2,1)  
mod.plot('D1', 'error', true)
hold on
mod.plot('equilibria');

% 
% xlim([L+100 R-100]);
% set(gca, 'Xtick', 0:1000:R-1);
% ylim([-30 10]);
% set(gca, 'Ytick', -30:15:10);

% text(gca, 0.92, 0.92, 'A', 'Units', 'normalized', 'fontsize', 15, 'tag','fignotext');
% set(findobj(gcf,'tag','fignotext'),'FontWeight','bold')
% set(gcf, 'position', [526.6000 400 400 300]);
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',18); % 刻度设置
xlabel('','fontsize',18,'fontname','Times')	%设置x轴字体
ylabel('Drift D^{(1)}','fontsize',18,'fontname','Times')	%设置y轴字体
set(gca, 'layer', 'top')
box on

% *************** *************** Diffusion: D_2 *************** ***************
subplot(3,2,2)  
mod.plot('D2', 'error', true)

% xlim([736 4227]);
% set(gca, 'Xtick', 0:1000:4000);
% ylim([-0.03 0.03]);
% set(gca, 'Ytick', -0.02:0.02:0.02);

% text(gca, 0.92, 0.92, 'B', 'Units', 'normalized', 'fontsize', 15, 'tag','fignotext');
% set(findobj(gcf,'tag','fignotext'),'FontWeight','bold')
% set(gcf, 'position', [526.6000 400 400 300]);
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',18); % 刻度设置
xlabel('','fontsize',18,'fontname','Times')	%设置x轴字体
ylabel('Diffusion D^{(2)}','fontsize',18,'fontname','Times')	%设置y轴字体
set(gca, 'layer', 'top')
box on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 稳态概率分布 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% P_st = A/D2(x)* exp{int(D1(x)/D2(x)}
subplot(3,2,3)
thepdf = mod.pdf;  % 考虑乘性噪声的影响，会导致吸引盆地的变化
mod.plot(thepdf);
% xpoint = linspace(0.28,0.75);
% pdff = thepdf.PDF_eq(xpoint);
% plot(xpoint,pdff)

% xlim([0.05, 0.55]);
% set(gca, 'Xtick', -5:5:5);
% ylim([0 10]);
% set(gca, 'Ytick', 0:2:10);
% text(gca, 0.92, 0.92, 'C', 'Units', 'normalized', 'fontsize', 15, 'tag','fignotext');
% set(findobj(gcf,'tag','fignotext'),'FontWeight','bold')
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
xlabel('','fontsize',18,'fontname','Times')	%设置x轴字体
ylabel('Probability distribution','fontsize',18,'fontname','Times')	%设置y轴字体
% set(gcf, 'position', [526.6000 400 400 300]);
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 有效势函数 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U_eff = ln(D2(x)) - int(D1(x)/D2(x)}
subplot(3,2,4)
mod.plot('potential_eff')  % 考虑乘性噪声的影响，会导致吸引盆地的变化

% xlim([0.1, 0.5]);
% set(gca, 'Xtick', -5:5:5);
% ylim([0 0.3]);
% set(gca, 'Ytick', 0:0.1:0.3);
% text(gca, 0.9, 0.92, 'D', 'Units', 'normalized', 'fontsize', 15, 'tag','fignotext');
% set(findobj(gcf,'tag','fignotext'),'FontWeight','bold')
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
xlabel('','fontsize',18,'fontname','Times')	%设置x轴字体
ylabel('Effective potential','fontsize',18,'fontname','Times')	%设置y轴字体
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 平均退出时间 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res = mod.mean_exit('all', thepdf);
fprintf('The weighted mean exit times:\n   the left basin : %g %s\n   the right basin : %g %s\n', res{1}.WT, mod.timeunit, res{2}.WT, mod.timeunit);

subplot(3,2,6)
mod.plot('mean_exit')
% 
% xlim([-5.5 5.5]);
% set(gca, 'Xtick', -5:5:5);
ylim([0 100]);
% set(gca, 'Ytick', 0:20:100);
% text(gca, 0.9, 0.92, 'F', 'Units', 'normalized', 'fontsize', 15, 'tag','fignotext');
% set(findobj(gcf,'tag','fignotext'),'FontWeight','bold') % 加粗
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
xlabel(mod.namex,'fontsize',18,'fontname','Times')	%设置x轴字体
ylabel('Mean exit time (month)','fontsize',18,'fontname','Times')	%设置y轴字体
% set(gcf, 'position', [526.6000 400 400 300]);
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 存活概率 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,5) 
mod.plot('survival','ylim',[0 150]);

% xlim([-5.5 5.5]);
% set(gca, 'Xtick', -5:5:5);
ylim([0 130]);
% set(gca, 'Ytick', 0:1000:3000);
% text(gca, 0.9, 0.92, 'E', 'Units', 'normalized', 'fontsize', 15, 'tag','fignotext');
% set(findobj(gcf,'tag','fignotext'),'FontWeight','bold') % 加粗
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
xlabel(mod.namex,'fontsize',18,'fontname','Times')	%设置x轴字体
ylabel('Survival time (month)','fontsize',18,'fontname','Times')	%设置y轴字体
set(gca, 'layer', 'top')
set(gcf, 'position', [300 -100 900 1200]);
box on

%%
close all;
% surv = mod.survival('all',[0 250]); 
subplot(2,1,1)
[~, res] = mod.plot('survival_func', mod.equilibria(3).x, surv);  
% 返回res中包括给定初值后的存活时间及对应存活概率

f = chebfun(res.surv, [res.t(1) res.t(end)], 'equi', 10000);  % 使用chebfun 插值得到更精确时间的存活概率
mediansurv = find(f == 0.5);  
fprintf('median exit time of initial Infectious I %f is %f\n',mod.equilibria(3).x, mediansurv)

plot([0 mediansurv mediansurv], [0.5 0.5 0], 'color', [137 23 23]/256, 'LineWidth', 1.8, 'LineStyle', '--');
xlim([0 250])

% plot([0 mediansurv mediansurv], [0.5 0.5 0], 'color', [137 23 23]/256, 'LineWidth', 1.8, 'LineStyle', '--');
% xlim([0 250])

set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
xlabel('Time','fontsize',18,'fontname','Times')	%设置x轴字体
ylabel('Probility of  survival time','fontsize',18,'fontname','Times')	%设置y轴字体
set(gca, 'layer', 'top')
% set(gcf, 'position', [526.6000 400 400 300]);
box on

%
subplot(2,1,2)
mod.plot('exit_distrib', mod.equilibria(3).x, surv);

% 计算平均退出时间
meansurv = cumsum(f);
meansurv = meansurv(meansurv.domain(2));
fprintf('mean exit time of initial Infectious I %f is %f\n', mod.equilibria(1).x, meansurv)

fS = -diff(f);  % 退出时间的概率密度函数 = -S'(xa,t)

hold on
plot([mediansurv mediansurv], [0 fS(mediansurv)],  'color', [137 23 23]/256, 'LineWidth', 1.8, 'LineStyle', '--'); % median exit time

hold on
plot([meansurv meansurv], [0 fS(meansurv)], 'k-', 'LineWidth', 1.8, 'LineStyle', '-');  % mean exit time
xlim([0 250])
% ylim([0 0.015])
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
xlabel('Time (month)','fontsize',18,'fontname','Times')	%设置x轴字体
ylabel('Density of exit times','fontsize',18,'fontname','Times')	%设置y轴字体
set(gca, 'layer', 'top')
set(gcf, 'position', [100 200 400 700]);
box on