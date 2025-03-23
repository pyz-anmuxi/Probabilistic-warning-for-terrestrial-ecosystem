
% Implements Langevin Reconstruction based on NDVI data
% 2025.2.21

clc;clear;close all;

%%
file = 'D:\CODE\Exit time\Ecosystem\Data_test\Data\';
ncFilePath = 'kNDVI_Example.nc';
var = ncread(strcat(file,ncFilePath),'var');

tab = var(5:end,2:2);
tim = 1: 1: size(tab);
xi = 1: 0.1: size(tab);
tab = interp1(tim, tab, xi, 'PCHIP'); % ���Բ�ֵ

%% Step2���ع�Langevin���̣��õ�D1��D2����
L = min(tab); 
R = max(tab); 

bins = 12;  % ����M_1�������������D_1���߹⻬�ȣ���ε��� 
% ���󣺻ع���ƾ����ȿ����޷��ﵽ��������
% ��С�����߲��⻬
Tau = 1:4;  % fitlm ���M(i)��Tau,Ӱ�����
% a list of time lags considered (integer values)
% (we chose 1:5. For high resolution data we can consider bigger number of time lags).
% x(t+tau) - x(t)
dt = 0.5; % time step of the data series Ӱ�������������
method = 'Nadaraya-Watson'; 
% you can choose if you want to estimate the cconditional moments 'directly' which is less smooth or via 'Nadaraya-Watson' estimator (our choice) which is smoother. 
% method = 'direct';

res = LangevinReconst(tab, L, R, bins, Tau, dt, method, 'kNDVI');

% save('ni','res', 'tim', 'tab', '-v7.3')

% Step3:ʹ�ò�ͬ�ķ��������ع��Ĳ�ȷ����
% % ����������ݿɲ���ִ��
% nreplicates = 10;  % set to zero for run without monte carlo error propagation
% 
% results = Langevin_bootstrap(data, res, {'error-propagation', 'block'}, nreplicates, struct('MaxStep', 0.01));
% % results.bootstrap������bootmethod��Ӧ��result��
% % reconstr���ع����ݵ㣩��result������õ���mean exit��pdf��
% save('data_bootstrap', 'results', 'tim', 'data' , '-v7.3')

% ��ȡD1 �� D2 ����
% nreplicates = 0; % set to zero for run without monte carlo error propagation
% if nreplicates > 0
%     results = load('Infectives_boot1');
% else
%     results = load('.\20230708\Measles.mat');
% end

% tim = results.tim;
% data = results.data;
% res = results.results;

mod = langevin_eq(res, 'weightedspline');

mod.namex = 'kNDVI';
mod.timeunit = '(year)';
% 
% fprintf('Small scale is %f\n', mod.equilibria(1).x)
% fprintf('scale threshold is %f\n', mod.equilibria(2).x)
% fprintf('big scale is %f\n', mod.equilibria(3).x * (max(tab_1) - min(tab_1))+ min(tab_1))
% *************** ***************Drift: D_1 *************** ***************
figure(2)  
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
set(gcf, 'position', [526.6000 400 400 300]);
set(gca,'linewidth',1.5)  % �߿�Ӵ�
set(gca,'FontName','Times New Roman','FontSize',18); % �̶�����
xlabel(mod.namex,'fontsize',18,'fontname','Times')	%����x������
ylabel('Drift D_1','fontsize',18,'fontname','Times')	%����y������
set(gca, 'layer', 'top')
box on

% *************** *************** Diffusion: D_2 *************** ***************
figure(3)   
mod.plot('D2', 'error', true)

% xlim([736 4227]);
% set(gca, 'Xtick', 0:1000:4000);
% ylim([-0.03 0.03]);
% set(gca, 'Ytick', -0.02:0.02:0.02);

% text(gca, 0.92, 0.92, 'B', 'Units', 'normalized', 'fontsize', 15, 'tag','fignotext');
% set(findobj(gcf,'tag','fignotext'),'FontWeight','bold')
set(gcf, 'position', [526.6000 400 400 300]);
set(gca,'linewidth',1.5)  % �߿�Ӵ�
set(gca,'FontName','Times New Roman','FontSize',18); % �̶�����
xlabel(mod.namex,'fontsize',18,'fontname','Times')	%����x������
ylabel('Diffusion D_2','fontsize',18,'fontname','Times')	%����y������
set(gca, 'layer', 'top')
box on

%% ʱ������ͼ

namex = 'kNDVI';
timeunit = '(year)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ʱ������ͼ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(xi, tab, 'Color',[130/255 177/255 255/255],'LineWidth',1.2)

% % hold on
line([1, 150], [mod.equilibria(1).x, mod.equilibria(1).x], 'LineWidth', 1.8, 'LineStyle', '--', 'color', [255/256 183/256 77/256]);
line([150, 200], [mod.equilibria(3).x, mod.equilibria(3).x], 'LineWidth', 1.8, 'LineStyle', '--', 'color',  [55/256 131/256 59/256]);
line([200, 396], [mod.equilibria(1).x, mod.equilibria(1).x],  'LineWidth', 1.8, 'LineStyle', '--', 'color', [255/256 183/256 77/256]);
% yline(mod.equilibria(1).x, '--',  'LineWidth', 1.8, 'color', [255/256 183/256 77/256]);
% yline(mod.equilibria(2).x,   'LineWidth', 1.8, 'color', [0.5 0.5 0.5]);
% yline(mod.equilibria(3).x,   'LineWidth', 1.8, 'color', [55/256 131/256 59/256]);

ylim([min(tab)-0.02 max(tab)+0.02]);
% set(gca, 'Ytick', 0:2500:5000)
set(gca, 'XTick',[3*12 13*12 23*12 33*12], 'XTickLabel', {'1985','1995','2005','2015'});
xlim([tim(1) tim(end)]);

% set(findobj(gcf, '-property', 'fontsize'), 'fontsize', 10);
set(gcf, 'position', [1200 500 600 240]);
% text(gca, 0.9, 0.95, 'A', 'Units', 'normalized', 'fontsize', 10, 'tag','fignotext');
set(gca,'linewidth',1.5)  % �߿�Ӵ�
set(gca,'FontName','Times New Roman','FontSize',15); % �̶�����
ylabel(namex,'fontsize',18,'fontname','Times')	%����x������
xlabel(['Time ' timeunit],'fontsize',18,'fontname','Times')	%����y������
box on

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ��̬���ʷֲ� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% P_st = A/D2(x)* exp{int(D1(x)/D2(x)}
figure(4)
thepdf = mod.pdf;  % ���ǳ���������Ӱ�죬�ᵼ��������صı仯
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
set(gca,'linewidth',1.5)  % �߿�Ӵ�
set(gca,'FontName','Times New Roman','FontSize',15); % �̶�����
xlabel(mod.namex,'fontsize',18,'fontname','Times')	%����x������
ylabel('Probability density','fontsize',18,'fontname','Times')	%����y������
set(gcf, 'position', [526.6000 400 400 300]);
box on

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ��Ч�ƺ��� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U_eff = ln(D2(x)) - int(D1(x)/D2(x)}
figure(5)  
mod.plot('potential_eff')  % ���ǳ���������Ӱ�죬�ᵼ��������صı仯

% xlim([0.1, 0.5]);
% set(gca, 'Xtick', -5:5:5);
% ylim([0 0.3]);
% set(gca, 'Ytick', 0:0.1:0.3);
% text(gca, 0.9, 0.92, 'D', 'Units', 'normalized', 'fontsize', 15, 'tag','fignotext');
% set(findobj(gcf,'tag','fignotext'),'FontWeight','bold')
set(gca,'linewidth',1.5)  % �߿�Ӵ�
set(gca,'FontName','Times New Roman','FontSize',15); % �̶�����
xlabel(mod.namex,'fontsize',18,'fontname','Times')	%����x������
ylabel('Effective potential','fontsize',18,'fontname','Times')	%����y������
set(gcf, 'position', [526.6000 400 400 300]);
box on

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ƽ���˳�ʱ�� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res = mod.mean_exit('all', thepdf);
fprintf('The weighted mean exit times:\n   the left basin : %g %s\n   the right basin : %g %s\n', res{1}.WT, mod.timeunit, res{2}.WT, mod.timeunit);

figure(6) 
mod.plot('mean_exit')
% 
% xlim([-5.5 5.5]);
% set(gca, 'Xtick', -5:5:5);
% ylim([0 100]);
% set(gca, 'Ytick', 0:20:100);
% text(gca, 0.9, 0.92, 'F', 'Units', 'normalized', 'fontsize', 15, 'tag','fignotext');
% set(findobj(gcf,'tag','fignotext'),'FontWeight','bold') % �Ӵ�
set(gca,'linewidth',1.5)  % �߿�Ӵ�
set(gca,'FontName','Times New Roman','FontSize',15); % �̶�����
xlabel(mod.namex,'fontsize',18,'fontname','Times')	%����x������
ylabel('Mean exit time','fontsize',18,'fontname','Times')	%����y������
set(gcf, 'position', [526.6000 400 400 300]);
box on

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ������ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7)  
mod.plot('survival','ylim',[0 300]);

% xlim([-5.5 5.5]);
% set(gca, 'Xtick', -5:5:5);
% ylim([0 3000]);
% set(gca, 'Ytick', 0:1000:3000);
% text(gca, 0.9, 0.92, 'E', 'Units', 'normalized', 'fontsize', 15, 'tag','fignotext');
% set(findobj(gcf,'tag','fignotext'),'FontWeight','bold') % �Ӵ�
set(gca,'linewidth',1.5)  % �߿�Ӵ�
set(gca,'FontName','Times New Roman','FontSize',15); % �̶�����
xlabel(mod.namex,'fontsize',18,'fontname','Times')	%����x������
ylabel('Survival time','fontsize',18,'fontname','Times')	%����y������
set(gca, 'layer', 'top')
set(gcf, 'position', [526.6000 400 400 300]);
box on

