%Plot the effect electrode-neuron distance for DBS for monopolar scenario 1
% i.e. the DBSMonopolar code (See details in that script)
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

ithresbemarray = [8.0687e-04 8.0371e-04 8.0617e-04 8.1108e-04 8.1363e-04 8.1600e-04 8.3074e-04...
                  8.5041e-04 8.5880e-04 8.7203e-04 9.3729e-04 0.0010 0.0011 0.0012...
                  0.0018 0.0030 0.0048 0.0074 0.0363 0.1223 0.3138 0.6587...
                  0.0101 0.0102 0.0102 0.0102 0.0102 0.0102 0.0102 0.0101 0.0101 0.0102...
                  0.0101 0.0101 0.0101 0.0100 0.0099 0.0108 0.0126 0.0154 0.0532 0.2240 0.6260 1.3595];
ithrescearray = [0.0011 0.0011 0.0011 0.0011 0.0011 0.0011 0.0011 0.0011 0.0011 0.0012 0.0013...
                 0.0015 0.0016 0.0018 0.0026 0.0041 0.0067 0.0073 0.0141 0.0438 0.1092 0.2269...
                 0.0041 0.0041 0.0041 0.0041 0.0041 0.0041 0.0041 0.0041 0.0041 0.0041 0.0042...
                 0.0042 0.0042 0.0043 0.0045 0.0050 0.0057 0.0066 0.0175 0.0698 0.1909 0.4127];
subplot(3,2,2)
semilogx(distance(1:14),ithresbemarray((1:14)+22),'LineStyle','-','Color',"#000080",'LineWidth', 3)% blue
hold on
semilogx(distance(1:14),ithrescearray((1:14)+22),'LineStyle','-.','Color',"#FFC000",'LineWidth', 3) %Gold
ax = gca();
ax.YColor = [0 0 0];
ylim([3*10^-3, 14*10^-3])
title("Monopolar DBS of a 4mm long, 2\mum diameter Axon", 'FontSize', 12);
error1_long = abs((-ithresbemarray((1:14)+22)+ithrescearray((1:14)+22))'./ithrescearray(23:44))*100;
legend('Bidomain BEM','Cable Equation','location','northwest')



subplot(3,2,4)
semilogx(distance(1:14),ithresbemarray((1:14)),'LineStyle','-','Color',"#000080",'LineWidth', 3)% blue
hold on
semilogx(distance(1:14),ithrescearray((1:14)),'LineStyle','-.','Color',"#FFC000",'LineWidth', 3) %Gold
ax = gca();
ax.YColor = [0 0 0];
ylim([7*10^-4, 18*10^-4])
error2_long = abs((-ithresbemarray(1:14)+ithrescearray(1:14))./ithrescearray(1:14))*100;
legend('Bidomain BEM','Cable Equation','location','northwest')

subplot(3,2,6)
semilogx(distance(1:14),error1_long(1:14),'LineStyle','-','Color',"#FF7518",'LineWidth', 3)% dark orange
hold on
semilogx(distance(1:14),error2_long(1:14),'LineStyle','-.','Color',"#FA8072",'LineWidth', 3) %light orange
ylim([0,220])
xlabel("d (m)", 'FontSize', 12);
legend("Perpendicular","Parallel",'location','northwest')

DBS.Units='inches';
DBS.Position=[-8.0104 4.9792 7 2];