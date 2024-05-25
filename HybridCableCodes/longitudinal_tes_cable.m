%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Boundary Element Method of Bidomain Modeling for Predicting Cellular
%                 Responses to Electromagnetic Fields
%
%             Code Authors: David Czerwonky, Luis Gomez
% Description: This script runs the hybrid cable approach to predict the
% behavior of HH axon embeded in a rectangular extracellular bath that is
% exposed to TES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('./Meshing/')
addpath('./Meshing/matrix')
addpath('./Meshing/mesh/')
addpath('./BEM_MEX_Fortran_code')
addpath('./channeldynamics/')
addpath('./Cableequation/')

amp = 11.2*10^-3;
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

%initalize time scale
tarr=0:dt:dur;
N=numel(tarr);
waveform=@(t) Iinject(t-.0002).*(t-.0002>0);
wave = waveform(tarr);

%define the locations of the axon compartments
p=zeros([M 3]);
p(:,3)=-len/2:len/(M-1):len/2;
p=p(:,[1 3 2]);
p(:,3)=p(:,3);
dl=p(2:end,:)-p(1:end-1,:);
peval=(p(2:end,:)+p(1:end-1,:))/2;


%calculate the right-hand side from the DBS electrode
rhs=TESRHS0(elecid,t2phead',phead');

%calculate the E-field produced by the DBS electrode
[Efield]=adjointBEM_electrode(t2phead',phead',epseff,rhs,peval');

%calculate the extracellular potentials
Ve=zeros([M 1]);
Ve(1)=0;
Ve(2:end)=cumsum(sum(dl'.*Efield,1));

%find the actication threshold (thval) and the transmembrane voltage (Vm)
Vm=multicompartsolv(a,len,rl,M,-amp,Ve,waveform,dt,dur);

%Spatial Plots of the transmembrane voltage
tarr = 1000*tarr;
Vm = 1000*Vm;
test_figure = figure;

%sample at 1/5th and 4/5th of the axon
id45 = ceil(4*M/5);
id15 = ceil(1*M/5);

xlabel("time (ms)",'FontSize', 12)
yyaxis left
hold on
plot(tarr,Vm(id15,:),'LineStyle','--','Color',"#ffda03 ",'LineWidth', 3) %Sunflower Yellow
plot(tarr,Vm(id45,:),'LineStyle','-.','Color'," #d6aa28",'LineWidth', 3)%Dijon
ylabel("Voltage (mV)",'FontSize', 12)
ylim([-80 50])
hold off
yyaxis right 
plot(tarr,1000*thval*wave,'-r','LineWidth', 3)
ylabel("Current (mA)",'FontSize', 12)
xlabel("time (ms)",'FontSize', 12)
test_figure.Units='inches';
test_figure.Position=[-8.0104 4.9792 2 2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                             End of Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%