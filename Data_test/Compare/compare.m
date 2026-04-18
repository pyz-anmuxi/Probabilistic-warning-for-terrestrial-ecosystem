clc; clear; close all;


file = 'D:\CODE\Exit time\Ecosystem\Data_test\Compare\';
data = load(strcat(file,'Mea_simuldata2.mat'));

xi = data.re_tab;
tim = 1:1:length(xi);

figure(1)
range = 1500:1:2000;
plot(tim(range), xi(range), 'Color',[130/255 177/255 255/255],'LineWidth',1.2)

set(gca, 'XTick',[range(150) range(300) range(450)], 'XTickLabel', {'150','300','450'});

set(gcf, 'position', [600 500 600 240]);
set(gca,'linewidth',1.5) 
set(gca,'FontName','Times New Roman','FontSize',15);  
ylabel('State V(t)','fontsize',18,'fontname','Times')	 
xlabel('The number of simulation data','fontsize',18,'fontname','Times')	 
box on

%%
clc; clear; close all;


file = 'D:\CODE\Exit time\Ecosystem\Data_test\Compare\';
data = load(strcat(file,'Mea_simuldata1.mat'));

xi = data.re_tab;
tim = 1:1:length(xi);

figure(3)
range = 5500:1:6000;
plot(tim(range), xi(range), 'Color',[130/255 177/255 255/255],'LineWidth',1.2)

set(gca, 'XTick',[range(150) range(300) range(450)], 'XTickLabel', {'150','300', '450'});

set(gcf, 'position', [600 500 600 240]);
set(gca,'linewidth',1.5)  
set(gca,'FontName','Times New Roman','FontSize',15); 
ylabel('State V(t)','fontsize',18,'fontname','Times')	 
xlabel('The number of simulation data','fontsize',18,'fontname','Times')	 
box on