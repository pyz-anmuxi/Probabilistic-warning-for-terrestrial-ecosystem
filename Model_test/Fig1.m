%% 根据模型实现Fig1中存活概率分布、平均退出时间
% 2025 
%% 
clc;clear;

addpath('D:\CODE\Exit time\chebfun-master');

asigma = 0.015;
% 设置真实模型为面向对象
realmod = openrealmod(asigma, [42 84]); % realmod: 1*1 langevin_eq

%% C图：存活概率S(x0,t)
% dS(x0,t)/dt = D1(x0) dS(x0,t)/dx0 + D2(x0) d^2S(x0,t)/dx0^2
% realmod.survival：pdepeq求解抛物型PDE
figure(1)

surv = realmod.survival('all',[0 3000]);  % 求解 S(x0,t) 时间范围（纵轴）[0 400] 取决于不同系统的不同存活时间
realmod.plot(surv, 'verticalbar', false);

hold on 
xline(realmod.equilibria(2).x, 'k--')
xline(realmod.equilibria(3).x, 'r--')

%% A图：平均退出时间 T(x0) 
% D1(x0) dT(x0,t)/dx0 + D2(x0) d^2T(x0,t)/dx0^2 = -1
% realmod.mean_exit：bvp45求解初边值问题
figure(2)

me = realmod.mean_exit('all', []);
realmod.plot(me, 'verticalbar', false);

hold on 
xline(realmod.equilibria(2).x, 'k--') % 不稳定的平衡点
xline(realmod.equilibria(3).x, 'b--')  % 第二个稳定平衡点
% ylim([0 3000])

%% B图：给定初始值的存活概率
% 找到给定初始值对应的存活概率
figure(3)

[~, res] = realmod.plot('survival_func', realmod.equilibria(3).x, surv);  
% 返回res中包括给定初值后的存活时间及对应存活概率

f = chebfun(res.surv, [res.t(1) res.t(end)], 'equi', 10000);  % 使用chebfun 插值得到更精确时间的存活概率
mediansurv = find(f == 0.5);  
fprintf('median exit time of initial Infectious I %f is %f\n', realmod.equilibria(3).x, mediansurv)

plot([0 mediansurv mediansurv], [0.5 0.5 0], 'r--');
xlim([0 4000])
set(gca,'layer','top')
box on
%% D图：给定初始值的退出时间
figure(4)

realmod.plot('exit_distrib', realmod.equilibria(3).x, surv);

% 计算平均退出时间
meansurv = cumsum(f);
meansurv = meansurv(meansurv.domain(2));
fprintf('mean exit time of initial Infectious I %f is %f\n', realmod.equilibria(3).x, meansurv)

fS = -diff(f);  % 退出时间的概率密度函数 = -S'(xa,t)

hold on
plot([mediansurv mediansurv], [0 fS(mediansurv)], 'r--');  % median exit time

hold on
plot([meansurv meansurv], [0 fS(meansurv)], 'b--'); % mean exit time
% xlim([0 400])
% ylim([0 0.015])
set(gca,'layer','top')
box on
