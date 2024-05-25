%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Boundary Element Method of Bidomain Modeling for Predicting Cellular
%                 Responses to Electromagnetic Fields
%
%             Code Author: David M. Czerwonky
%
% Description: This script runs the hybrid cable approach to predict the
% behavior of HH axon embeded in an infinite extracellular space that is
% exposed to DBS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load supporting software
addpath('./Meshing/')

%load axon properties
a=10^-6;   %axon radius
len=.004;  %axon length
rl=1;      %axon longitudinal resistivity
M=100;     %number of compartments

%time scale
dt=1*10^-6;     %time step
dur=20*10^-3; % duration

%define electrode waveform
Iinject=@(t) (t/(5*10^-5).*(t<5*10^-5)+(t>=5*10^-5)).*(t<5*10^-3);

%define DBS electrode geometry
Na=8;
Nl=100;
[p_dbs,tr_dbs]=cylmesh(Na,Nl,.01,.001/2,0);
p_dbs=p_dbs(:,[1 3 2]);
tr_dbs=tr_dbs(:,[1,3,2]);
londir=2;
p_dbs(:,londir)=p_dbs(:,londir)-.003;
p_dbs(:,1)=p_dbs(:,1)-(.001/2)-(10^-6)-2.5*10^-6;
pcen_dbs=(p_dbs(tr_dbs(:,1),:)+p_dbs(tr_dbs(:,2),:)+p_dbs(tr_dbs(:,3),:))/3;
sigmaext=2; %extracelluar conductivity
elecid_dbs=zeros(size(pcen_dbs(:,3)));
elecid_dbs((pcen_dbs(:,londir)>.00075).*(pcen_dbs(:,londir)<.00225)==1)=1/sigmaext;
elecid_dbs((pcen_dbs(:,londir)<-.00075).*(pcen_dbs(:,londir)>-.00225)==1)=-1/sigmaext;

%define epsilon effective for finding the E-fields later on
epseff=(0-sigmaext)/(0+sigmaext)*ones([numel(tr_dbs)/3 1]);%(-sigma_in+sigma_out)/(sigma_in-sigma_out)

%find the actication threshold (thval) and the transmembrane voltage (Vm)
[thval Vm waveform tarr]=dbs_head(a,len,rl,Iinject,p_dbs',tr_dbs',elecid_dbs,epseff,dt,dur,M);

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
plot(tarr,1000*thval*wave,'-r','LineWidth', 3)
ylabel("Current (mA)",'FontSize', 10)
xlabel("time (ms)",'FontSize', 10)
test_figure.Units='inches';
test_figure.Position=[-8.0104 4.9792 2 2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                             End of Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%