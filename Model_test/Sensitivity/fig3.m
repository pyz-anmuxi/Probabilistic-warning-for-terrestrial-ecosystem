%
% figure 3 Sensitivity analysis of data size and some parameters: dt, ndata, bins and tau
%

clc;clear;close all;
asigma = 0.019;
realmod = openrealmod(asigma, [41 85]);

%% here a 10x10 analysis of the effect of dt and ndata
% takes about 24 hours to run
dts = 0.1:0.1:1.0;  % 保存数据的时间精度
ndata = 60000:10000:150000;  % 数据量
nreplicates = 50;

%10x10 I can improve that later
results = cell(length(dts), length(ndata));

for i = 1:length(dts)
    for j = 1:length(ndata)
        fprintf('Running dt = %g, ndata = %g\n', dts(i), ndata(j));
        results{i, j} = Langevin_bootstrap([], realmod, 'error-propagation' , nreplicates, ...
        struct('MaxStep', 0.005, 'RemoveBias', true, 'datasize', [1, ndata(j)], 'domain', ...
        realmod.domain, 'nbins', 50, 'Tau', 1:5, 'dt', dts(i), 'method', 'Nadaraya-Watson', 'title', 'Figure x - data resolution'));
    end
    save('Fig_resol', 'results', 'dts', 'ndata');
end

%% 读取数据绘图
file = '.\';
results_resol = load(strcat(file,'Fig_resol'));
dts = results_resol.dts;
ndata = results_resol.ndata;
results_resol = results_resol.results;
        
[ndatas, dts1] = meshgrid(ndata, dts);
ranges_av1 = zeros(size(dts1));
ranges_av2 = zeros(size(dts1));
biasD1 = zeros(size(dts1));
biasD2 = zeros(size(dts1));
        
for i = 1:size(dts1, 1)
    for j = 1:size(dts1, 2)
        res = getbootresults(results_resol{i, j});
        ranges_av1(i, j) = res.relrange(1);
        ranges_av2(i, j) = res.relrange(2);
        biasD1(i, j) = res.relbias(1);
        biasD2(i, j) = res.relbias(2);
    end
end

%% %%%%%%%%%%%%%% % dt, ndata 对重构误差的敏感性分析
figure(1)

C = [068 004 090
     065 062 133
     048 104 141
     031 146 139
     145 213 066
     248 230 032] / 255;
C = flip(C, 1); % 反转矩阵的行顺序

n = 256; % 生成 256 种颜色
C_interp = interp1(1:size(C, 1), C, linspace(1, size(C, 1), n));
colormap(C_interp);

X1_normalized = (biasD1 - min(biasD1(:))) / (max(biasD1(:)) - min(biasD1(:)));
% colormap( [linspace(1, 0.5, 50)', linspace(1, 0.5, 50)', ones(50, 1)]);
[~,~] = contourf(dts1, ndatas / 1000, X1_normalized, 8);
hold on
plot(0.2, 80, 'k*');

% xlim([0.1 1.0])
% set(gca, 'xtick', 0.1:0.2:1.0)
set(gca, 'ytick', 60:30:150)

title('Bias D_1')
xlabel('Data resolution ({\Delta}t)','fontsize',15,'fontname','Times');
ylabel('Number of points (x1000)','fontsize',15,'fontname','Times');

set(gcf, 'position', [600 500 410 450]);
set(gca,'linewidth',1.2)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
% text(gca, 0.9, 0.95, 'G', 'Units', 'normalized', 'fontsize', 10);
shading flat
c = colorbar;
% set(c,'YTick',0.012:0.004:0.024); %色标值范围及显示间隔

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
colormap(C_interp);
X2_normalized = (biasD2 - min(biasD2(:))) / (max(biasD2(:)) - min(biasD2(:)));
[~,~] = contourf(dts1, ndatas / 1000, X2_normalized, 8);
hold on
plot(0.2, 80, 'k*');
% 
% xlim([0.1 1.0])
% set(gca, 'xtick', 0.1:0.2:1.0)
set(gca, 'ytick', 60:30:150)

title('Bias D_2')
xlabel('Data resolution ({\Delta}t)','fontsize',15,'fontname','Times');
ylabel('Number of points (x1000)','fontsize',15,'fontname','Times');

set(gcf, 'position', [600 500 410 450]);
set(gca,'linewidth',1.2)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
% text(gca, 0.9, 0.95, 'G', 'Units', 'normalized', 'fontsize', 10);
shading flat
c = colorbar;
% set(c,'YTick',0.004:0.004:0.016); %色标值范围及显示间隔

%% %%%%%%%%%%%%%% % dt, ndata 对平均退出时间的不确定性分析
figure(3)
ranges_av1(ranges_av1 > 1) = 1;
% min(min(ranges_av1))  max(max(ranges_av1))
colormap(C_interp);
contourf(dts1, ndatas/1000, ranges_av1, 8);
hold on
plot(0.2, 80, 'k*');

% xlim([0.1 1.0])
% set(gca, 'xtick', 0.1:0.2:1.0)
set(gca, 'ytick', 60:30:150)

% text(gca, 0.9, 0.95, 'E', 'Units', 'normalized', 'fontsize', 10);
title('Uncertainty (domain 1)')
xlabel('Data resolution ({\Delta}t)','fontsize',15,'fontname','Times');
ylabel('Number of points (x1000)','fontsize',15,'fontname','Times');

set(gcf, 'position', [600 500 410 450]);
set(gca,'linewidth',1.2)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
shading flat
c = colorbar;
% set(c,'YTick',0.1:0.2:1); %色标值范围及显示间隔

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
ranges_av2(ranges_av2 > 1) = 1;
% min(min(ranges_av2))  max(max(ranges_av2))
colormap(C_interp);
contourf(dts1, ndatas/1000, ranges_av2, 8);
hold on
plot(0.2, 80, 'k*');
% 
% xlim([0.1 1.0])
% set(gca, 'xtick', 0.1:0.2:1.0)
set(gca, 'ytick', 60:30:150)

% text(gca, 0.9, 0.95, 'E', 'Units', 'normalized', 'fontsize', 10);
title('Uncertainty (domain 2)')
xlabel('Data resolution ({\Delta}t)','fontsize',15,'fontname','Times');
ylabel('Number of points (x1000)','fontsize',15,'fontname','Times');

set(gcf, 'position', [600 500 410 450]);
set(gca,'linewidth',1.2)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
shading flat
c = colorbar;
set(c,'YTick',0.3:0.1:1); %色标值范围及显示间隔

%% here a 5x5 analysis of the effect of bins and tau
% 重构D1和D2时，bin和tau的选择
% takes about 6 hours to run
nbinss = 10:10:70;  % 重构时的点数
maxtau = 3:2:12;  % 
nreplicates = 50;

Ndata = 80000;  % 生成时间序列的数据数量
dt = 0.2;  % 时间精度

% create nreplicates*2条 data
res = realmod.simulate(linspace(0, Ndata * dt, Ndata), realmod.draw_pdf(nreplicates * 2, 1), struct('MaxStep', 0.005));
save('res','res');

res = load( 'res.mat');
res = res.res';
x = res.y;

results = cell(length(nbinss), length(maxtau));  % 创建5*5的矩阵
for i = 1:length(nbinss)
    for j = 1:length(maxtau)
        fprintf('Running nbins = %g, maxtau = %g\n', nbinss(i), maxtau(j));
        results{i, j} = Langevin_bootstrap([], realmod, 'replicates' , nreplicates, ...
        struct('MaxStep', 0.005, 'data', x, 'RemoveBias', true, 'datasize', size(x), 'domain', ...
        realmod.domain, 'bins', nbinss(i), 'Tau', 1:maxtau(j), 'dt', dt, 'method', 'Nadaraya-Watson', 'title', 'Figure x - nbins'));
    end
    
end
save('Fig_bins', 'results', 'nbinss', 'maxtau', 'x');

%% 读取数据绘图
file = '.\';
results = load(strcat(file,'Fig_bins'));
nbinss = results.nbinss;
maxtau = results.maxtau;
results_bin = results.results;

[maxtaus, nbinss1] = meshgrid(maxtau, nbinss);
ranges_av1 = zeros(size(nbinss1));  % 区域1加权平均退出时间的不确定性
ranges_av2 = zeros(size(nbinss1));   % 区域2加权平均退出时间的不确定性
biasD1 = zeros(size(nbinss1));  % D1 的bias
biasD2 = zeros(size(nbinss1));

for i = 1:size(nbinss1, 1)
    for j = 1:size(nbinss1, 2)
        res = getbootresults(results_bin{i, j});
        ranges_av1(i, j) = res.relrange(1);
        ranges_av2(i, j) = res.relrange(2);
        biasD1(i, j) = res.relbias(1);
        if length(res.relbias) > 1
            biasD2(i, j) = res.relbias(2);
        end
    end
end
       
%% %%%%%%%%%%%%%% % bin 和tau 对重构误差的敏感性分析%% 
figure(5)


C = [068 004 090
     065 062 133
     048 104 141
     031 146 139
     145 213 066
     248 230 032] / 255;
C = flip(C, 1); % 反转矩阵的行顺序

n = 256; % 生成 256 种颜色
C_interp = interp1(1:size(C, 1), C, linspace(1, size(C, 1), n));
colormap(C_interp);

X1_normalized = (biasD1 - min(biasD1(:))) / (max(biasD1(:)) - min(biasD1(:)));
contourf(maxtaus, nbinss1, X1_normalized, 8)
xlim([min(maxtau) max(maxtau)])
hold on
plot(5, 50, 'k*');

xlim([min(maxtau) max(maxtau)])
set(gca, 'xtick', 3:2:12)
set(gca, 'ytick', 10:20:70)

% text(gca, 0.9, 0.95, 'E', 'Units', 'normalized', 'fontsize', 10);
title('BiasD_1')
xlabel('Max \tau','fontsize',15,'fontname','Times');
ylabel('Number of bins','fontsize',15,'fontname','Times');

set(gcf, 'position', [600 500 410 450]);
set(gca,'linewidth',1.2)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
shading flat
c = colorbar;
% set(c,'YTick',0.013:0.0015:0.022); %色标值范围及显示间隔

%%%%%%%%%%%%%%%%%
figure(6) 
colormap(C_interp);

X2_normalized = (biasD2 - min(biasD2(:))) / (max(biasD2(:)) - min(biasD2(:)));
contourf(maxtaus, nbinss1, X2_normalized, 8)
xlim([min(maxtau) max(maxtau)])
hold on
plot(5, 50, 'k*');

xlim([min(maxtau) max(maxtau)])
set(gca, 'xtick', 3:2:12)
set(gca, 'ytick', 10:20:70)

% text(gca, 0.9, 0.95, 'E', 'Units', 'normalized', 'fontsize', 10);
title('BiasD_2')
xlabel('Max \tau','fontsize',15,'fontname','Times');
ylabel('Number of bins','fontsize',15,'fontname','Times');

set(gcf, 'position', [600 500 410 450]);
set(gca,'linewidth',1.2)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
shading flat
c = colorbar;
% set(c,'YTick',0.1:0.2:1); %色标值范围及显示间隔

%% %%%%%%%%%%%%%% % bin 和tau 对退出时间的敏感性分析%%  
figure(7)

colormap(C_interp);

rang1_normalized = (ranges_av1 - min(ranges_av1(:))) / (max(ranges_av1(:)) - min(ranges_av1(:)));
contourf(maxtaus, nbinss1, rang1_normalized, 8)
hold on
plot(5, 50, 'k*');

xlim([min(maxtau) max(maxtau)])
set(gca, 'xtick', 3:2:12)
set(gca, 'ytick', 10:20:70)

% text(gca, 0.9, 0.95, 'E', 'Units', 'normalized', 'fontsize', 10);
title('Uncertainty (domain 2)')
xlabel('Max \tau','fontsize',15,'fontname','Times');
ylabel('Number of bins','fontsize',15,'fontname','Times');

set(gcf, 'position', [600 500 410 450]);
set(gca,'linewidth',1.2)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
shading flat
c = colorbar;
% set(c,'YTick',0.1:0.2:1); %色标值范围及显示间隔

%%%%%%%%%%%%%%%%%%
figure(8)

colormap(C_interp);

rang2_normalized = (ranges_av2 - min(ranges_av2(:))) / (max(ranges_av2(:)) - min(ranges_av2(:)));
contourf(maxtaus, nbinss1, rang2_normalized, 8)
hold on
plot(5, 50, 'k*');

xlim([min(maxtau) max(maxtau)])
set(gca, 'xtick', 3:2:12)
set(gca, 'ytick', 10:20:70)

% text(gca, 0.9, 0.95, 'E', 'Units', 'normalized', 'fontsize', 10);
title('Uncertainty (domain 1)')
xlabel('Max \tau','fontsize',15,'fontname','Times');
ylabel('Number of bins','fontsize',15,'fontname','Times');

set(gcf, 'position', [600 500 410 450]);
set(gca,'linewidth',1.2)  % 边框加粗
set(gca,'FontName','Times New Roman','FontSize',15); % 刻度设置
shading flat
c = colorbar;
% set(c,'YTick',1.2:0.01:1.25); %色标值范围及显示间隔

