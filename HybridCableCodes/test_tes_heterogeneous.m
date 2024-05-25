%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Boundary Element Method of Bidomain Modeling for Predicting Cellular
%                 Responses to Electromagnetic Fields
%
%             Code Author: David M. Czerwonky
%
% Description: This script runs the hybrid cable approach for the
% HH axon embeded in the heterogeneous extracellular bath of section 3.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('./BEM_MEX_Fortran_code')
addpath('./channeldynamics/')
addpath('./Cableequation/')
addpath('./Meshing/')
addpath('./Meshing/mesh/')
addpath('./Meshing/matrix/')

rl=1;
sigmaext=0.6;

%Threshold guess
ampval = 10^-3;

%electrode waveform information
close all
Iinject=@(t) (t/(5*10^-5).*(t<5*10^-5)+(t>=5*10^-5)).*(t<5*10^-3);
dur1=10*10^-2;
dur=20*10^-3;
dt=1*10^-5;
tarr=0:dt:dur;
N=numel(tarr);
wave=zeros(size(tarr));
wave(1:numel(tarr))=Iinject(tarr-0.0002).*(tarr>0.0002);

%% Initalize Geometry
load('./Meshing/2LayerRectBath_res10.mat','p','tri','id','V','F','V2','F2','V3','F3')
% trisurf(tri,p(:,1),p(:,2),p(:,3),'edgealpha',0.2,'facealpha',0.2,'facecolor','#FF0000')
p_sur = p;
clear p

%% Define material parameters
sig=[1/rl 0.276 1.654];
sigav=zeros(numel(id),1);
sigdiff=zeros(numel(id),1);

%define sigma average
sigav(id==1) = (sig(1)+sig(2))/2;
sigav(id==2) = (sig(2)+sig(3))/2;
sigav(id==3) = (sig(3)+0)/2;

%define sigma difference
sigdiff(id==1) = (sig(2)-sig(1));
sigdiff(id==2) = (sig(3)-sig(2));
sigdiff(id==3) = (0-sig(3));

%define epsilon effective
m = numel(F)/3;
epseff = -sigdiff((m+1):end)./(2*sigav((m+1):end));

%create heterogenerous extracellular bath
[phead,t2phead] = cat_meshes(V2,F2,V3,F3);

%find compartment locations
M = 100;
len = 4*10^-3;
a = 10^-6;
p=zeros([M 3]);
p(:,3)=0:len/(M-1):len;
p=p(:,[1 3 2]);
dl=p(2:end,:)-p(1:end-1,:);
peval=(p(2:end,:)+p(1:end-1,:))/2;

%create electrode vector
elecid2 = zeros([numel(tri)/3 1]);
%find the center of each triangle element
pcen=(phead(t2phead(:,1),:)+phead(t2phead(:,2),:)+phead(t2phead(:,3),:))/3;
elecdir=2; % set the electrodes along the y-axis [1 is x-axis 3 is z-axis]
for i=1:numel(pcen)/3
    %1.654 is the extracellular conductivity
    if(pcen(i,elecdir)<=min(pcen(:,elecdir))+.000001)
        elecid2(i) = 1/1.654; %find the contribution from the cathode
    elseif(pcen(i,elecdir)>=max(pcen(:,elecdir))-.000001)
        elecid2(i) = -1/1.654; %find the contribution from the cathode
    end
end
%create the forcing term/ right-hand side of the matrix eqn
rhs=DBSRHS0(elecid2,t2phead',phead');

%Find the Extracellular Potentials from the E-field
[Efield]=runcodet2prhs(t2phead',phead',epseff,rhs,peval');
Ve=zeros([M 1]); %initalize extracellular potentials
Ve(1)=0;
Ve(2:end)=cumsum(sum(dl'.*Efield,1));

%set helper functions for bisection algorithm
getRow = @(data, rowNum) data(rowNum,10:N); % Helper function
thfunc=@(amp) max(getRow(multicompartsolv(a,len,rl,M,-amp,Ve,Iinject,dt,dur),round(M/5)))>0;

%Use bisection method to find the activation threshold
thval=findthreshold(@(amp) thfunc(amp),1,10^-3,.1)

%find the transmembrane voltage waveforms at the activation threshold
Vm=multicompartsolv(a,len,rl,M,-thval,Ve,Iinject,dt,dur);
max(max(Vm))

%Spatial Plots of the transmembrane voltage
tarr = 1000*tarr;
Vm = 1000*Vm;
test_figure = figure;

%sample at 1/5th and 4/5th of the axon
id31 = ceil(4*M/5);
id32 = ceil(1*M/5);

xlabel("time (ms)",'FontSize', 10)
yyaxis left
hold on
plot(tarr,Vm(id31,:),'LineStyle','-','Color',"#FFC000",'LineWidth', 3) %Gold
plot(tarr,Vm(id32,:),'LineStyle','-.','Color',"#FFC000",'LineWidth', 3) %Gold
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