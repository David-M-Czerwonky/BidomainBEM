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
subplot(3,4,1)
loglog(distance,ithresbemarray,'LineStyle','-','Color',"#000080",'LineWidth', 3)% blue
hold on
load('Transverse_DBSelecposCE.mat')
loglog(distance,ithrescearray','LineStyle','-.','Color',"#FFC000",'LineWidth', 3) %Gold
ax = gca();
ax.YColor = [0 0 0];
ylim([10^-3,10^3])
% xlabel("d (m)", 'FontSize', 10);
ylabel("Perpendicular Threshold (A)", 'FontSize', 12);
title("DBS of a 4mm long, 2\mum diameter Axon", 'FontSize', 12);
error1_std = abs((-ithresbemarray+ithrescearray)./ithrescearray)*100;
legend('Bidomain BEM','Cable Equation','location','northwest')

load('Longitudinal_DBSelecposBEM.mat')
subplot(3,4,5)
loglog(distance,ithresbemarray,'LineStyle','-','Color',"#000080",'LineWidth', 3)% blue
hold on
load('Longitudinal_DBSelecposCE.mat')
loglog(distance,ithrescearray,'LineStyle','-.','Color',"#FFC000",'LineWidth', 3) %Gold
ax = gca();
ax.YColor = [0 0 0];
ylim([10^-4,10^0])
% xlabel("d (m)", 'FontSize', 10);
ylabel("Parallel Threshold (A)", 'FontSize', 11);
% title("Parallel Orientation 4mm long, 2\mum diameter", 'FontSize', 10);
legend('Bidomain BEM','Cable Equation','location','northwest')


subplot(3,4,9)
error2_std = abs((-ithresbemarray+ithrescearray)./ithrescearray)*100;
semilogx(distance,error1_std,'LineStyle','-','Color',"#FF7518",'LineWidth', 3)% dark orange
hold on
semilogx(distance,error2_std,'LineStyle','-.','Color',"#FA8072",'LineWidth', 3) %light orange
ylim([0,100])
xlabel("d (m)", 'FontSize', 12);
ylabel("Percentatge Difference (%)", 'FontSize', 12);
% title("% Difference 4mm long, 2\mum diameter", 'FontSize', 10);
legend("Perpendicular","Parallel",'location','northwest')

min(error1_std)
min(error2_std)

%% Long Axon Size

load('Transverse_DBSLong.mat')
subplot(3,4,2)
loglog(distance,ithresbemarray,'LineStyle','-','Color',"#000080",'LineWidth', 3)% blue
hold on
loglog(distance,ithrescearray','LineStyle','-.','Color',"#FFC000",'LineWidth', 3) %Gold
ax = gca();
ax.YColor = [0 0 0];
ylim([10^-3,10^3])
% xlabel("d (m)", 'FontSize', 10);
% ylabel("Activation Threshold (A)", 'FontSize', 10);
title("DBS of a 10mm long, 2\mum diameter Axon", 'FontSize', 12);
error1_long = abs((-ithresbemarray+ithrescearray)./ithrescearray)*100;
legend('Bidomain BEM','Cable Equation','location','northwest')



load('Longitudinal_DBSLong.mat')
subplot(3,4,6)
loglog(distance,ithresbemarray,'LineStyle','-','Color',"#000080",'LineWidth', 3)% blue
hold on
loglog(distance,ithrescearray(1:22),'LineStyle','-.','Color',"#FFC000",'LineWidth', 3) %Gold
ax = gca();
ax.YColor = [0 0 0];
ylim([10^-4,10^0])
% xlabel("d (m)", 'FontSize', 10);
% ylabel("Activation Threshold (A)", 'FontSize', 10);
% title("Parallel Orientation 10mm long, 2\mum diameter", 'FontSize', 10);
error2_long = abs((-ithresbemarray+ithrescearray(1:22))./ithrescearray(1:22))*100;
legend('Bidomain BEM','Cable Equation','location','northwest')

min(error1_long)
min(error2_long)

subplot(3,4,10)
semilogx(distance,error1_long,'LineStyle','-','Color',"#FF7518",'LineWidth', 3)% dark orange
hold on
semilogx(distance,error2_long,'LineStyle','-.','Color',"#FA8072",'LineWidth', 3) %light orange
ylim([0,100])
xlabel("d (m)", 'FontSize', 12);
% ylabel("Percentatge (%)", 'FontSize', 10);
% title("% Difference 10mm long, 2\mum diameter", 'FontSize', 10);
legend("Perpendicular","Parallel",'location','northwest')

DBS.Units='inches';
DBS.Position=[-8.0104 4.9792 7 2];