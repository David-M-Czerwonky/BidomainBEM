%Plot the effect electrode-neuron distance for DBS
addpath('E:\MATLAB\NeuronBEM\Results')
addpath('E:\MATLAB\NeuronBEM\CableBEM\')
addpath('C:\Users\david\Documents\MATLAB Bulky Files')
close all

distance = [50*10^-9, 80*10^-9, 10^-7, 5*10^-7, 8*10^-7, 10^-6, 2.5*10^-6, 5*10^-6, 7.5*10^-6,...
	 10^-5, 2.5*10^-5, 5*10^-5, 7.5*10^-5, 10^-4, 2.5*10^-4, 5*10^-4, 7.5*10^-4, 10^-3,...
     2.5*10^-3, 5*10^-3, 7.5*10^-3, 10^-2];

%% Standard Axon Size
DBS = figure();
load('Transverse_DBSelecposBEM.mat')
subplot(3,3,1)
semilogx(distance(1:14),ithresbemarray(1:14),'LineStyle','-','Color',"#000080",'LineWidth', 3)% blue
hold on
load('Transverse_DBSelecposCE.mat')
semilogx(distance(1:14),ithrescearray(1:14)','LineStyle','-.','Color',"#FFC000",'LineWidth', 3) %Gold
ax = gca();
ax.YColor = [0 0 0];
ylim([0.14,0.22])
% xlabel("d (m)", 'FontSize', 10);
ylabel("Perpendicular Threshold (A)", 'FontSize', 12);
title("4mm long, 2\mum diameter", 'FontSize', 12);
error1_std = abs((-ithresbemarray+ithrescearray)./ithrescearray)*100;
legend('Bidomain BEM','Cable Equation','location','northwest')

load('Longitudinal_DBSelecposBEM.mat')
subplot(3,3,4)
loglog(distance(1:14),ithresbemarray(1:14),'LineStyle','-','Color',"#000080",'LineWidth', 3)% blue
hold on
load('Longitudinal_DBSelecposCE.mat')
semilogx(distance(1:14),ithrescearray(1:14),'LineStyle','-.','Color',"#FFC000",'LineWidth', 3) %Gold
ax = gca();
ax.YColor = [0 0 0];
% ylim([10^-4,10^0])
% xlabel("d (m)", 'FontSize', 10);
ylabel("Parallel Threshold (A)", 'FontSize', 11);
% title("Parallel Orientation 4mm long, 2\mum diameter", 'FontSize', 10);
legend('Bidomain BEM','Cable Equation','location','northwest')


subplot(3,3,7)
error2_std = abs((-ithresbemarray+ithrescearray)./ithrescearray)*100;
semilogx(distance(1:14),error1_std(1:14),'LineStyle','-','Color',"#FF7518",'LineWidth', 3)% dark orange
hold on
semilogx(distance(1:14),error2_std(1:14),'LineStyle','-.','Color',"#FA8072",'LineWidth', 3) %light orange
ylim([0,33])
xlabel("d (m)", 'FontSize', 12);
ylabel("Error (%)", 'FontSize', 12);
% title("% Difference 4mm long, 2\mum diameter", 'FontSize', 10);
legend("Perpendicular","Parallel",'location','northwest')

min(error1_std)
min(error2_std)

%% Long Axon Size

load('Transverse_DBSLong.mat')
subplot(3,3,2)
semilogx(distance(1:14),ithresbemarray(1:14),'LineStyle','-','Color',"#000080",'LineWidth', 3)% blue
hold on
semilogx(distance(1:14),ithrescearray(1:14)','LineStyle','-.','Color',"#FFC000",'LineWidth', 3) %Gold
ax = gca();
ax.YColor = [0 0 0];
ylim([0.16,0.23])
% xlabel("d (m)", 'FontSize', 10);
% ylabel("Activation Threshold (A)", 'FontSize', 10);
title("10mm long, 2\mum diameter", 'FontSize', 12);
error1_long = abs((-ithresbemarray+ithrescearray)./ithrescearray)*100;
legend('Bidomain BEM','Cable Equation','location','northwest')



load('Longitudinal_DBSLong.mat')
subplot(3,3,5)
semilogx(distance(1:14),ithresbemarray(1:14),'LineStyle','-','Color',"#000080",'LineWidth', 3)% blue
hold on
semilogx(distance(1:14),ithrescearray(1:14),'LineStyle','-.','Color',"#FFC000",'LineWidth', 3) %Gold
ax = gca();
ax.YColor = [0 0 0];
% ylim([10^-4,10^0])
% xlabel("d (m)", 'FontSize', 10);
% ylabel("Activation Threshold (A)", 'FontSize', 10);
% title("Parallel Orientation 10mm long, 2\mum diameter", 'FontSize', 10);
error2_long = abs((-ithresbemarray+ithrescearray(1:22))./ithrescearray(1:22))*100;
legend('Bidomain BEM','Cable Equation','location','northwest')

min(error1_long)
min(error2_long)

subplot(3,3,8)
semilogx(distance(1:14),error1_long(1:14),'LineStyle','-','Color',"#FF7518",'LineWidth', 3)% dark orange
hold on
semilogx(distance(1:14),error2_long(1:14),'LineStyle','-.','Color',"#FA8072",'LineWidth', 3) %light orange
ylim([0,20])
xlabel("d (m)", 'FontSize', 12);
% ylabel("Percentatge (%)", 'FontSize', 10);
% title("% Difference 10mm long, 2\mum diameter", 'FontSize', 10);
legend("Perpendicular","Parallel",'location','northwest')


%% Thick Axon Size

load('Transverse_DBSThick5um.mat')
subplot(3,3,3)
semilogx(distance(1:14),ithresbemarray(1:14),'LineStyle','-','Color',"#000080",'LineWidth', 3)% blue
hold on
semilogx(distance(1:14),ithrescearray(1:14)','LineStyle','-.','Color',"#FFC000",'LineWidth', 3) %Gold
ax = gca();
ax.YColor = [0 0 0];
ylim([0.03,0.1])
% xlabel("d (m)", 'FontSize', 10);
% ylabel("Activation Threshold (A)", 'FontSize', 10);
title("4mm long, 10\mum diameter", 'FontSize', 12);
error1_thick = abs((-ithresbemarray+ithrescearray)./ithrescearray)*100;
legend('Bidomain BEM','Cable Equation','location','northwest')



load('Longitudinal_DBSThick5um.mat')
subplot(3,3,6)

semilogx(distance(1:14),ithresbemarray(1:14),'LineStyle','-','Color',"#000080",'LineWidth', 3)% blue
hold on
semilogx(distance(1:14),ithrescearray(1:14),'LineStyle','-.','Color',"#FFC000",'LineWidth', 3) %Gold
ax = gca();
ax.YColor = [0 0 0];
% ylim([10^-4,10^0])
% xlabel("d (m)", 'FontSize', 10);
% ylabel("Activation Threshold (A)", 'FontSize', 10);
% title("Parallel Orientation 4mm long, 10\mum diameter", 'FontSize', 10);
error2_thick = abs((-ithresbemarray+ithrescearray(1:22))./ithrescearray(1:22))*100;
legend('Bidomain BEM','Cable Equation','location','northwest')

min(error1_thick)
min(error2_thick)

subplot(3,3,9)
semilogx(distance(1:14),error1_thick(1:14),'LineStyle','-','Color',"#FF7518",'LineWidth', 3)% dark orange
hold on
semilogx(distance(1:14),error2_thick(1:14),'LineStyle','-.','Color',"#FA8072",'LineWidth', 3) %light orange
ylim([0,73])
xlabel("d (m)", 'FontSize', 12);
% ylabel("Percentatge (%)", 'FontSize', 10);
% title("% Difference 4mm long, 10\mum diameter", 'FontSize', 10);
legend("Perpendicular","Parallel",'location','northwest')

DBS.Units='inches';
DBS.Position=[-8.0104 4.9792 7 2];