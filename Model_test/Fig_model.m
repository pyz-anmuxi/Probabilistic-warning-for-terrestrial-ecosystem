%% vegetation model generated data
% 
% 2023.9.14

clc;clear;close all;

%% 设置面向对象

%addpath('D:\Exit time\code\matlab-code\chebfun-master');

realmod = [];
asigma = 0.019;

% openrealmod：设置面向对象obj = langevin_eq的properties：D1,D2,DD2,equilibria,domain,namex
if isempty(realmod)
    realmod = openrealmod(asigma, [41 85]);  % realmod:
end

realmod.plot(realmod.pdf)   % 判断能否产生状态切换数据
%% Step1: create data from the model
% 若已经生成数据，可以不用执行！！
% simulate(obj, tspan,y0,options)使用Euler-Maruyama方法对langevin方程时间积分
                           %--obj = 1*1 langevin_eq (面向对象)
                           %--tspan：积分的时间跨度
                           %y0 = drawa-pdf(1) (累计密度函数在插值点处的值)
                           %--obj.odefun--创建ODE函数：D1(x) + sqrt(2*dt*D2(x))*randn(1)/dt
                           %--simulate_rng--euler
                                   
% can be set to true to generate new data
% needs grind for matlab to be installed: download from https://www.sparcs-center.org/grind

% stabilizing is not needed if we draw from the stationary pdf
Ndata = 80000; % 数据量
dt = 0.2; % 保存数据的时间精度
% MaxStep：积分精度

% 模型生成数据
res = realmod.simulate(linspace(0, Ndata * dt, Ndata), realmod.draw_pdf(1), struct('MaxStep', 0.005)); %struct: 使用名为字段的数据容器将相关数据组合在一起的数据类型。

x = res.y;
t = res.t;

data = table(t, x);
writetable(data, 'model_data.csv')
% 查看生成数据
plot(t,x,'b-')
xlabel('Time');
ylabel('Tree cover (%)');
%% Step2：重构Langevin方程，得到D1和D2
file = '.\';
tab = readtable(strcat(file,'model_data.csv'));

L = realmod.domain(1); % min(tab.x)
R = realmod.domain(2); % max(tab.x)
bin = 50;
Tau = 1:5;

res = LangevinReconst(tab.x, L, R, bin, Tau, tab.t(2) - tab.t(1), 'Nadaraya-Watson', 'modelled data');

% save('model_Reconst', 'res' , '-v7.3');

%% 比较原模型、一组数据下重构数据点、Langevin重构后的 D1 和 D2 
mod = langevin_eq(res, 'weightedspline');
% d1 = fit(res.C, res.D1);
% D1 = chebfun(d1(l), [L1 R1], 'equi', 5000);
x = res.C;

D1 = res.D1;
D2 = res.D2;

ax1 = subplot(2,1,1);
hax = realmod.plot('D1', 'verticalbar', false, 'error', false);  % 实际方程 realmod
set(hax, 'Color', 'r', 'LineWidth', 1, 'LineStyle', ':');
hold on
hax2 = mod.plot('D1');  % 重构得到Langevin方程的结果
hax1 = plot(x, D1, 'b--','linewidth',1.5);  % 重构后的数据结果
ylabel('Drift D_1')
legend([hax1,hax2,hax],'Reconst','Langevin\_Reconst','Real\_model');

ax2 = subplot(2,1,2);
h3 = realmod.plot('D2', 'verticalbar', false, 'error', false);
set(h3, 'Color', 'r', 'LineWidth', 1, 'LineStyle', ':');
hold on
h2 = mod.plot('D2');
h1 = plot(x,D2, 'b--','linewidth',1.5);
ylabel('Drift D_2');
legend([h1, h2, h3],'Reconst','Langevin\_Reconst','Real\_model');

xlabel('State x')
xticklabels(ax1,{}) % move the xticklabels of subfigure(1)

%% Step3: 使用不同的方法评估重构的不确定性
% % 若已经生成数据，可以不用执行！！
nreplicates = 50; % 复制次数，与MC方法有关

results = Langevin_bootstrap(tab.x, res, {'model-based'}, nreplicates, struct('MaxStep', 0.005, 'model', realmod));
% tab.x：模型生成的时间序列数据
% res：step2中的重构的 D1，D2 等数据
% model-based：基于原始模型生成bootstrap data
% error-propagation:基于重构方程生成数据集
% block：块自举法，简单百分位复制法

save('model_bootstrap', 'results' , '-v7.3');
% model_bootstrap.mat相比于model_Reconst.mat多结构体bootstrap  
% 结构体bootstrap包括每组生成数据的reconstr(重构结果)和result(计算得到的pdfs和mean exit time)
%% Step4: 计算平稳分布，有效势，退出时间，存活时间

nreplicates = 10;
if nreplicates > 0
    results = load(strcat(file,'model_bootstrap'));
else
    results = load(strcat(file,'model_Reconst'));
end

results = results. results;
mod = langevin_eq(results, 'weightedspline');
mod.namex = 'Tree cover (%)';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 时间序列图 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tab = readtable(strcat(file,'model_data.csv'));
Ndata = 80000; % 数据量

figure(10) 
timextick = 0:5000:Ndata;
plot(tab.t, tab.x, 'Color',  [196 165 222]./256)
% [196 165 222]
% 
% xlabel(['Time (day)')
% ylabel(mod.namex)
% ylim(mod.domain);  % 只保留大于零的
set(gca, 'Xtick', timextick);
% xlim([min(tab.t), max(tab.t)]);

set(gcf, 'position', [526.6000 631.4000 800 300]);
% text(gca, 0.92, 0.92, 'A', 'Units', 'normalized', 'fontsize', 15, 'tag','fignotext');
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',18); % 刻度设置
xlabel('Time (day)', 'fontsize',18,'fontname','Times')	%设置x轴字体
ylabel('Tree cover (%)','fontsize',18,'fontname','Times')	%设置y轴字体
set(gca,'Layer','top')
box on

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Drift: D_1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)  
mod.plot('D1', 'error', true)
hold on
mod.plot('equilibria');

hold on
hax = realmod.plot('D1', 'verticalbar', false, 'error', false);  % 真实模型
set(hax, 'Color', [47/256 31/256 73/256], 'LineWidth', 1.2, 'LineStyle', '--');

xlim(mod.domain);
ylim([-0.6 0.4])

set(gcf, 'position', [526.6000 400 400 300]);
% text(gca, 0.92, 0.92, 'B', 'Units', 'normalized', 'fontsize', 15, 'tag','fignotext');
% set(findobj(gcf,'tag','fignotext'),'FontWeight','bold')
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',18); % 刻度设置
xlabel(mod.namex,'fontsize',18,'fontname','Times')	%设置x轴字体
ylabel('Drift D_1','fontsize',18,'fontname','Times')	%设置y轴字体
set(gca,'Layer','top')
box on

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Diffusion: D_2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)   
mod.plot('D2', 'error', true)

hold on
hax = realmod.plot('D2', 'verticalbar', false, 'error', false);
set(hax, 'Color', [47/256 31/256 73/256], 'LineWidth', 1.2, 'LineStyle', '--');

xlim(mod.domain);
ylim([0.2 1.2])
% set(gca, 'Ytick', 0.2:0.015:0.045);

set(gcf, 'position', [526.6000 400 400 300]);
% text(gca, 0.92, 0.92, 'C', 'Units', 'normalized', 'fontsize', 15, 'tag','fignotext');
% set(findobj(gcf,'tag','fignotext'),'FontWeight','bold')
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',18); % 刻度设置
xlabel(mod.namex,'fontsize',18,'fontname','Times')	%设置x轴字体
ylabel('Diffusion D_2','fontsize',18,'fontname','Times')	%设置y轴字体
set(gca,'Layer','top')
box on
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 稳态概率分布 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure(3)
thepdf = mod.pdf;
mod.plot(thepdf);

% ylim([0 0.1])
% set(gca, 'Ytick', 0:0.2:1);

set(gcf, 'position', [526.6000 400 400 300]);
% text(gca, 0.92, 0.92, 'D', 'Units', 'normalized', 'fontsize', 15, 'tag','fignotext');
% set(findobj(gcf,'tag','fignotext'),'FontWeight','bold')
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
xlabel(mod.namex,'fontsize',18,'fontname','Times')	%设置x轴字体
ylabel('Probability distribution','fontsize',18,'fontname','Times')	%设置y轴字体
box on

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 有效势函数 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)  
mod.plot('potential_eff')

set(gcf, 'position', [526.6000 400 400 300]);
% text(gca, 0.92, 0.92, 'E', 'Units', 'normalized', 'fontsize', 15, 'tag','fignotext');
% set(findobj(gcf,'tag','fignotext'),'FontWeight','bold')
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',18); % 刻度设置
xlabel(mod.namex,'fontsize',18,'fontname','Times')	%设置x轴字体
ylabel('Effective potential ','fontsize',18,'fontname','Times')	%设置y轴字体
box on
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 平均退出时间 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res = mod.mean_exit('all', thepdf);
fprintf('The weighted mean exit times:\n   the left basin: %g %s\n   the right basin : %g %s\n', res{1}.WT, mod.timeunit, res{2}.WT, mod.timeunit);

figure(5) 
mod.plot('mean_exit')

res1 = realmod.mean_exit('all');
hax = realmod.plot(res1, 'verticalbar', false, 'error', true);
set(hax, 'Color', [255 236 76]/256, 'LineWidth', 1.8, 'LineStyle', '--');

ylim([0 800])
set(gca, 'Ytick', 0:200:800);

set(gcf, 'position', [526.6000 400 400 300]);
% text(gca, 0.92, 0.92, 'F', 'Units', 'normalized', 'fontsize', 15, 'tag','fignotext');
% set(findobj(gcf,'tag','fignotext'),'FontWeight','bold')
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',18); % 刻度设置
xlabel(mod.namex,'fontsize',18,'fontname','Times')	%设置x轴字体
ylabel('Mean exit time (day)','fontsize',18,'fontname','Times')	%设置y轴字体
set(gca,'Layer','top')
box on

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 存活概率 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6)  
mod.plot('survival',[0 1000]);
set(gca, 'Ytick', 0:300:1000);
set(gcf, 'position', [526.6000 400 400 300]);
% text(gca, 0.92, 0.92, 'G', 'Units', 'normalized', 'fontsize', 15, 'tag','fignotext');
% set(findobj(gcf,'tag','fignotext'),'FontWeight','bold')
set(gca,'linewidth',1.5)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',18); % 刻度设置
xlabel(mod.namex,'fontsize',18,'fontname','Times')	%设置x轴字体
ylabel('Survival time (day)','fontsize',18,'fontname','Times')	%设置y轴字体
set(gca,'Layer','top')
box on
