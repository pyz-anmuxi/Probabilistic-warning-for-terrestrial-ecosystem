
% Implements Langevin Reconstruction based on k-NDVI data
% 2025.4.15

clc;clear;close all;

%%
file = 'D:\CODE\Exit time\Ecosystem\Data_test\Data\';
ncFilePath = 'kNDVI_Example.nc';
var = ncread(strcat(file,ncFilePath),'var');

tab = var(5:end,2:2); % Continuous sudden shift

% Data preprocessing: Z-score normalization
tab = (tab - mean(tab)) / std(tab);

tim = 1: 1: size(tab);
xi = 1: 0.5: size(tab);
tab = interp1(tim, tab, xi, 'PCHIP');  % Linear interpolation

%% Step2£ºReconstruct the Langevin equation to obtain D1 and D2 data.
L = min(tab); 
R = max(tab); 

bins = 9;  % Determine the amount of data in M_1, the smoothness of the fitted D_1 curve

Tau = 1:5;  % fitlm fits M(i) and Tau, affecting the error.
% a list of time lags considered (integer values)
% (we chose 1:5. For high resolution data we can consider bigger number of time lags).
% x(t+tau) - x(t)
dt = 0.5; % time step of the data series 
method = 'Nadaraya-Watson'; 
% you can choose if you want to estimate the cconditional moments 'directly' which is less smooth or via 'Nadaraya-Watson' estimator (our choice) which is smoother. 
% method = 'direct';

res = LangevinReconst(tab, L, R, bins, Tau, dt, method, 'kNDVI');

% save('ni','res', 'tim', 'tab', '-v7.3')

% check

mod = langevin_eq(res, 'weightedspline');

mod.namex = 'kNDVI';
mod.timeunit = '(year)';
mod.plot('D1', 'error', true)
hold on
mod.plot('equilibria');
plot(res.C,res.D1)

figure(3)   
mod.plot('D2', 'error', true)


%% Time series with metastable states

namex = 'Vegetation state, V';
timeunit = '(month)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Time series %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(xi, tab, 'Color',[130/255 177/255 255/255],'LineWidth',1.2)

% % hold on
hig = mod.equilibria(3).x;
low = mod.equilibria(1).x;
line([1, 150], [low low], 'LineWidth', 1.8, 'LineStyle', '--', 'color', [0.6350, 0.0780, 0.1840]);
line([150, 200], [hig hig], 'LineWidth', 1.8, 'LineStyle', '--', 'color',  [55/256 131/256 59/256]);
line([200, 396], [low low],  'LineWidth', 1.8, 'LineStyle', '--', 'color', [0.6350, 0.0780, 0.1840]);

ylim([min(tab)-0.02 max(tab)+0.02]);
% set(gca, 'Ytick', 0:2500:5000)
set(gca, 'XTick',[3*12 13*12 23*12 33*12], 'XTickLabel', {'1985','1995','2005','2015'});
xlim([tim(1) tim(end)]);

set(gcf, 'position', [1200 500 600 240]);
set(gca,'linewidth',1.5)  
set(gca,'FontName','Times New Roman','FontSize',15); 
ylabel(namex,'fontsize',18,'fontname','Times')	
xlabel(['Time ' timeunit],'fontsize',18,'fontname','Times')	
box on

%% data from D1 and D2

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


set(gca,'linewidth',1.5)  
set(gca,'FontName','Times New Roman','FontSize',18);
xlabel('','fontsize',18,'fontname','Times')
ylabel('Drift D^{(1)}','fontsize',18,'fontname','Times')
set(gca, 'layer', 'top')
box on

% *************** *************** Diffusion: D_2 *************** ***************
subplot(3,2,2)  
mod.plot('D2', 'error', true)

set(gca,'linewidth',1.5) 
set(gca,'FontName','Times New Roman','FontSize',18);
xlabel('','fontsize',18,'fontname','Times')	
ylabel('Diffusion D^{(2)}','fontsize',18,'fontname','Times')	
set(gca, 'layer', 'top')
box on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% probability distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% P_st = A/D2(x)* exp{int(D1(x)/D2(x)}
subplot(3,2,3)
thepdf = mod.pdf;
mod.plot(thepdf);
% xpoint = linspace(0.28,0.75);
% pdff = thepdf.PDF_eq(xpoint);
% plot(xpoint,pdff)


set(gca,'linewidth',1.5)  
set(gca,'FontName','Times New Roman','FontSize',15); 
xlabel('','fontsize',18,'fontname','Times')	
ylabel('Probability distribution','fontsize',18,'fontname','Times')	
% set(gcf, 'position', [526.6000 400 400 300]);
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% potential landscape %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U_eff = ln(D2(x)) - int(D1(x)/D2(x)}
subplot(3,2,4)
mod.plot('potential_eff')  


set(gca,'linewidth',1.5)  
set(gca,'FontName','Times New Roman','FontSize',15); 
xlabel('','fontsize',18,'fontname','Times')	
ylabel('Effective potential','fontsize',18,'fontname','Times')	
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mean exit time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res = mod.mean_exit('all', thepdf);
fprintf('The weighted mean exit times:\n   the left basin : %g %s\n   the right basin : %g %s\n', res{1}.WT, mod.timeunit, res{2}.WT, mod.timeunit);

subplot(3,2,6)
mod.plot('mean_exit')
% 
% xlim([-5.5 5.5]);
% set(gca, 'Xtick', -5:5:5);
ylim([0 100]);

set(gca,'linewidth',1.5)  
set(gca,'FontName','Times New Roman','FontSize',15);
xlabel(mod.namex,'fontsize',18,'fontname','Times')	
ylabel('Mean exit time (month)','fontsize',18,'fontname','Times')
% set(gcf, 'position', [526.6000 400 400 300]);
box on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% survival probability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,2,5) 
mod.plot('survival','ylim',[0 150]);

% xlim([-5.5 5.5]);
% set(gca, 'Xtick', -5:5:5);
ylim([0 130]);
% set(gca, 'Ytick', 0:1000:3000);

set(gca,'linewidth',1.5) 
set(gca,'FontName','Times New Roman','FontSize',15); 
xlabel(mod.namex,'fontsize',18,'fontname','Times')
ylabel('Survival time (month)','fontsize',18,'fontname','Times')
set(gca, 'layer', 'top')
set(gcf, 'position', [300 -100 900 1200]);
box on

%%
close all;
% surv = mod.survival('all',[0 250]); 
subplot(2,1,1)
[~, res] = mod.plot('survival_func', mod.equilibria(3).x, surv);  


f = chebfun(res.surv, [res.t(1) res.t(end)], 'equi', 10000); 
mediansurv = find(f == 0.5);  
fprintf('median exit time of initial Infectious I %f is %f\n',mod.equilibria(3).x, mediansurv)

plot([0 mediansurv mediansurv], [0.5 0.5 0], 'color', [137 23 23]/256, 'LineWidth', 1.8, 'LineStyle', '--');
xlim([0 250])

% plot([0 mediansurv mediansurv], [0.5 0.5 0], 'color', [137 23 23]/256, 'LineWidth', 1.8, 'LineStyle', '--');
% xlim([0 250])

set(gca,'linewidth',1.5)  
set(gca,'FontName','Times New Roman','FontSize',15);
xlabel('Time','fontsize',18,'fontname','Times')	
ylabel('Probility of  survival time','fontsize',18,'fontname','Times')
set(gca, 'layer', 'top')
% set(gcf, 'position', [526.6000 400 400 300]);
box on

%
subplot(2,1,2)
mod.plot('exit_distrib', mod.equilibria(3).x, surv);


meansurv = cumsum(f);
meansurv = meansurv(meansurv.domain(2));
fprintf('mean exit time of initial Infectious I %f is %f\n', mod.equilibria(1).x, meansurv)

fS = -diff(f);  % -S'(xa,t)

hold on
plot([mediansurv mediansurv], [0 fS(mediansurv)],  'color', [137 23 23]/256, 'LineWidth', 1.8, 'LineStyle', '--'); % median exit time

hold on
plot([meansurv meansurv], [0 fS(meansurv)], 'k-', 'LineWidth', 1.8, 'LineStyle', '-');  % mean exit time
xlim([0 250])
% ylim([0 0.015])
set(gca,'linewidth',1.5) 
set(gca,'FontName','Times New Roman','FontSize',15); 
xlabel('Time (month)','fontsize',18,'fontname','Times')	
ylabel('Density of exit times','fontsize',18,'fontname','Times')
set(gca, 'layer', 'top')
set(gcf, 'position', [100 200 400 700]);
box on