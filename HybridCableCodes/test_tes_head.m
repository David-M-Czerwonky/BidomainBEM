%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Boundary Element Method of Bidomain Modeling for Predicting Cellular
%                 Responses to Electromagnetic Fields
%
%             Code Authors: David Czerwonky, Luis Gomez
% Description: This script runs the hybrid cable approach to predict the
% behavior of HH axon embeded in a rectangular extracellular bath that is
% exposed to TES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load supporting software
addpath('./Meshing/')

%load axon properties
a=1*10^-6;   %axon radius
len=.004;  %axon length
rl=1;      %axon longitudinal resistivity
M=100;     %number of compartments

%time scale
dt=4*10^-6;     %time step
dur=20*10^-3; % duration

%define electrode waveform
Iinject=@(t) (t/(5*10^-5).*(t<5*10^-5)+(t>=5*10^-5)).*(t<5*10^-3);

%define extracellular bath geometry
[phead,t2phead,Q] = cube(10,10,10);
phead(:,1)=phead(:,1)-1/2;
phead(:,2)=phead(:,2)-1/2;
phead(:,3)=phead(:,3)-1/2;
phead(:,1)=.016*phead(:,1);
phead(:,2)=.008*phead(:,2);
phead(:,3)=.016*phead(:,3);

%define epsilon effective for finding the E-fields later on
sigmaext=2;
epseff=(sigmaext-0)/(sigmaext+0)*ones([numel(t2phead)/3 1]);%(sigma_in-sigma_out)/(sigma_in-sigma_out)

%define electrode positions
pcen=(phead(t2phead(:,1),:)+phead(t2phead(:,2),:)+phead(t2phead(:,3),:))/3;
elecid=zeros([numel(pcen)/3 1]);
elecdir=2; % 1 and 3 are transverse, 2 is longitudinal
elecid(pcen(:,elecdir)>=max(pcen(:,elecdir)))=1/sigmaext;
elecid(pcen(:,elecdir)<=min(pcen(:,elecdir)))=-1/sigmaext;

[thval Vm wave tarr]=tes_head(a,len,rl,Iinject,phead',t2phead',elecid,epseff,dt,dur,M);

%Spatial Plots of the transmembrane voltage
tarr = 1000*tarr;
Vm = 1000*Vm;
test_figure = figure;

%sample at 1/5th and 4/5th of the axon
id45 = ceil(4*M/5);
id15 = ceil(1*M/5);

xlabel("time (ms)",'FontSize', 10)
yyaxis left
hold on
plot(tarr,Vm(id15,:),'LineStyle','--','Color',"#ffda03 ",'LineWidth', 3) %Sunflower Yellow
plot(tarr,Vm(id45,:),'LineStyle','-.','Color'," #d6aa28",'LineWidth', 3)%Dijon
ylabel("Voltage (mV)",'FontSize', 10)
ylim([-80 50])
hold off
yyaxis right 
plot(tarr,1000*thval*wave,'-r','LineWidth', 3)
ylabel("Current (mA)",'FontSize', 10)
xlabel("time (ms)",'FontSize', 10)
test_figure.Units='inches';
test_figure.Position=[-8.0104 4.9792 2 2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                             End of Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%