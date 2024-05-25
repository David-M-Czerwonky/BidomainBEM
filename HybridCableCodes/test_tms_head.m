%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Boundary Element Method of Bidomain Modeling for Predicting Cellular
%                 Responses to Electromagnetic Fields
%
%             Code Author: David M. Czerwonky
%
% Description: This script runs the hybrid cable approach to predict the
% behavior of HH axon embeded in a spherical extracellular bath that is
% exposed to TMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load supporting software
addpath('./Meshing/')

%generated extracellular space
sigmaext=2;
[phead,t2phead] = icosphere(3);
phead=phead*.085;

%define epsilon effective for finding the E-fields later on
epseff=(sigmaext-0)/(sigmaext+0)*ones([numel(t2phead)/3 1]); %(sigma_in-sigma_out)/(sigma_in+sigma_out)

%load coil source information
load coilinformation.mat;
rv(:,3)=rv(:,3)+.09; 

%load axon properties
a=10^-6;   %axon radius
len=.004;  %axon length
rl=1;      %axon longitudinal resistivity
M=100;     %number of compartments

%time scale
dt=10^-6;     %time step
dur=20*10^-3; % duration

%find the activation threshold and transmembrane voltage Vm
[thval Vm waveform tarr]=tms_head(a,len,rl,VL,L,Jdc',rv',phead',t2phead',epseff,dt,dur,M);

%Spatial Plots of the transmembrane voltage
wave = waveform(tarr);
tarr = 1000*tarr;
Vm = 1000*Vm;
test_figure = figure;

%sample at 1/5th and 4/5th of the axon
id15 = ceil(1*M/5);
id45 = ceil(4*M/5);

%begin plotting
xlabel("time (ms)",'FontSize', 10)
yyaxis left
hold on
plot(tarr,Vm(id15,:),'LineStyle','-','Color',"#FFC000",'LineWidth', 3) %Gold
plot(tarr,Vm(id45,:),'LineStyle','-.','Color',"#FFC000",'LineWidth', 3) %Gold
ylabel("Voltage (mV)",'FontSize', 10)
ylim([-80 50])
hold off
yyaxis right 
plot(tarr,thval*wave/(10^9),'-r','LineWidth', 3)
ylabel("Current (kA/\mu s)",'FontSize', 10)
xlabel("time (ms)",'FontSize', 10)
test_figure.Units='inches';
test_figure.Position=[-8.0104 4.9792 2 2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                             End of Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%