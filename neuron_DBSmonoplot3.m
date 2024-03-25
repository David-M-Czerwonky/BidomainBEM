%Plot the effect electrode-neuron distance for DBS for monopolar scenario 1
% i.e. the DBSMonopolar3 code (See details in that script)
close all

distance = [50*10^-9, 80*10^-9, 10^-7, 5*10^-7, 8*10^-7, 10^-6, 2.5*10^-6, 5*10^-6, 7.5*10^-6,...
	 10^-5, 2.5*10^-5, 5*10^-5, 7.5*10^-5, 10^-4, 2.5*10^-4, 5*10^-4, 7.5*10^-4, 10^-3,...
     2.5*10^-3, 5*10^-3, 7.5*10^-3, 10^-2];

%% Standard Axon Size
DBS = figure();
load('Transverse_DBSelecposBEM.mat')
subplot(3,2,1)
semilogx(distance(1:14),ithresbemarray(1:14),'LineStyle','-','Color',"#000080",'LineWidth', 3)% blue
hold on
load('Transverse_DBSelecposCE.mat')
semilogx(distance(1:14),ithrescearray(1:14)','LineStyle','-.','Color',"#FFC000",'LineWidth', 3) %Gold
ax = gca();
ax.YColor = [0 0 0];
ylim([0.14 0.22])
% xlabel("d (m)", 'FontSize', 10);
ylabel("Perpendicular Threshold (A)", 'FontSize', 12);
title("Bipolar DBS of a 4mm long, 2\mum diameter Axon", 'FontSize', 12);
error1_std = abs((-ithresbemarray+ithrescearray)./ithrescearray)*100;
legend('Bidomain BEM','Cable Equation','location','northwest')

load('Longitudinal_DBSelecposBEM.mat')
subplot(3,2,3)
semilogx(distance(1:14),ithresbemarray(1:14),'LineStyle','-','Color',"#000080",'LineWidth', 3)% blue
hold on
load('Longitudinal_DBSelecposCE.mat')
loglog(distance(1:14),ithrescearray(1:14),'LineStyle','-.','Color',"#FFC000",'LineWidth', 3) %Gold
ax = gca();
ax.YColor = [0 0 0];
ylim([7*10^-4, 18*10^-4])
ylabel("Parallel Threshold (A)", 'FontSize', 11);
legend('Bidomain BEM','Cable Equation','location','northwest')


subplot(3,2,5)
error2_std = abs((-ithresbemarray+ithrescearray)./ithrescearray)*100;
semilogx(distance(1:14),error1_std(1:14),'LineStyle','-','Color',"#FF7518",'LineWidth', 3)% dark orange
hold on
semilogx(distance(1:14),error2_std(1:14),'LineStyle','-.','Color',"#FA8072",'LineWidth', 3) %light orange
ylim([0,30])
xlabel("d (m)", 'FontSize', 12);
ylabel("Percentatge Difference (%)", 'FontSize', 12);
legend("Perpendicular","Parallel",'location','northwest')

%% Monopolar Simulation

ithresbemarray = [0.0019 0.0019 0.0019 0.0019 0.0020 0.0020 0.0020...
                  0.0021 0.0021 0.0021 0.0023 0.0026 0.0028 0.0031...
                  0.0041 0.0061 0.0085 0.0115 0.0458 0.2097 0.6031 1.3313...
                  0.0011 0.0011 0.0011 0.0011 0.0011 0.0011 0.0011 0.0011 0.0011 0.0012...
                  0.0012 0.0013 0.0014 0.0015 0.0021 0.0034 0.0052 0.0076 0.0418 0.2063 0.5974 1.3219];
ithrescearray = [4.1981e-04 4.1981e-04 4.2349e-04 4.2349e-04 4.2126e-04 4.2262e-04 4.2641e-04...
                 4.2864e-04 4.3294e-04 4.3548e-04 4.5007e-04...
                 4.8147e-04 5.0623e-04 5.3053e-04 7.0188e-04...
                 9.9811e-04 9.9811e-04 9.9811e-04 9.9811e-04 9.9550e-04 9.9511e-04 0.0010...
                 9.9955e-04 0.0010 0.0010 0.0010 0.0010 0.0011 0.0012 0.0014 0.0017 0.0021 0.0042...
                 0.0042 0.0042 0.0043 0.0045 0.0050 0.0057 0.0066 0.0072 0.0322 0.1909 0.4127];
subplot(3,2,2)
semilogx(distance(1:14),ithresbemarray((1:14)+22),'LineStyle','-','Color',"#000080",'LineWidth', 3)% blue
hold on
semilogx(distance(1:14),ithrescearray((1:14)+22),'LineStyle','-.','Color',"#FFC000",'LineWidth', 3) %Gold
ax = gca();
ax.YColor = [0 0 0];
ylim([0.5*10^-3, 5*10^-3])
title("Monopolar DBS of a 4mm long, 2\mum diameter Axon", 'FontSize', 12);
error1_long = abs((-ithresbemarray((1:14)+22)+ithrescearray((1:14)+22))'./ithrescearray(23:44))*100;
legend('Bidomain BEM','Cable Equation','location','northwest')



subplot(3,2,4)
semilogx(distance(1:14),ithresbemarray((1:14)),'LineStyle','-','Color',"#000080",'LineWidth', 3)% blue
hold on
semilogx(distance(1:14),ithrescearray((1:14)),'LineStyle','-.','Color',"#FFC000",'LineWidth', 3) %Gold
ax = gca();
ax.YColor = [0 0 0];
ylim([1*10^-4, 4*10^-3])
error2_long = abs((-ithresbemarray(1:14)+ithrescearray(1:14))./ithrescearray(1:14))*100;
legend('Bidomain BEM','Cable Equation','location','northwest')

subplot(3,2,6)
semilogx(distance(1:14),error1_long(1:14),'LineStyle','-','Color',"#FF7518",'LineWidth', 3)% dark orange
hold on
semilogx(distance(1:14),error2_long(1:14),'LineStyle','-.','Color',"#FA8072",'LineWidth', 3) %light orange
xlabel("d (m)", 'FontSize', 12);
legend("Perpendicular","Parallel",'location','northwest')

DBS.Units='inches';
DBS.Position=[-8.0104 4.9792 7 2];