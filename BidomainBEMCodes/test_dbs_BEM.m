%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Boundary Element Method of Bidomain Modeling for Predicting Cellular
%                 Responses to Electromagnetic Fields
%
%             Code Author: David M. Czerwonky
%
% Description: This script finds the threshold of a DBS of an HH axon 4mm
% away and parallel to the electrode.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Supporting Software
addpath('./FortranBEMFix/matlab')
addpath('./FortranBEMFix/matlab/NeuronScripts/')
addpath('./channeldynamics/')
addpath('./Meshing/')

close all
clear all
clc

%% Initalize the time scale
dt=4*10^-6;     %time step
dur=20*10^-3; % duration
tarr=0:dt:dur;
N=numel(tarr);

%% Initalize the DBS electrode model
%generate electrode waveform
Iinject=@(t) (t/(5*10^-5).*(t<5*10^-5)+(t>=5*10^-5)).*(t<5*10^-3); %  5ms square wave pulse
wave=zeros(size(tarr)); %waveform
wave(1:numel(tarr))=Iinject(tarr-0.0002).*(tarr>0.0002); %delay pulse by 0.2 ms

%generate a simple bipolar DBS electrode
th = pi/2; % parallel orientation
% th = 0; % perpendicular orientation
dist = 4*10^-3; %distance between the neuron and DBS electrode on the x-axis
Na=8;
Nl=100;
[p_dbs,tr_dbs]=cylmesh(Na,Nl,.01,.001/2,0);
londir=3;
p_dbs(:,londir)=p_dbs(:,londir)-.003;
p_dbs(:,1)=p_dbs(:,1)-.001/2-10^-6-dist;
pcen_dbs=(p_dbs(tr_dbs(:,1),:)+p_dbs(tr_dbs(:,2),:)+p_dbs(tr_dbs(:,3),:))/3;
%define cathode and annode locations in the electrode array
elecid_dbs=zeros(size(pcen_dbs(:,3)));
elecid_dbs((pcen_dbs(:,londir)>.00075).*(pcen_dbs(:,londir)<.00225)==1)=1;
elecid_dbs((pcen_dbs(:,londir)<-.00075).*(pcen_dbs(:,londir)>-.00225)==1)=-1;
%rotate the DBS electrode
ptemp(:,1)=p_dbs(:,2)*cos(th)-p_dbs(:,3)*sin(th);
ptemp(:,2)=p_dbs(:,2)*sin(th)+p_dbs(:,3)*cos(th);
p_dbs(:,2)=ptemp(:,1);
p_dbs(:,3)=ptemp(:,2);
ndbs=numel(tr_dbs)/3;

%% Initalize Geometry
% generate neuron geometry
Na=8;
Nl=100;
[p_in,tr_in]=cylmesh(Na,Nl,.004,10^-6,0);
id1=20*Na;
id2=40*Na;
id3=60*Na;
id4=80*Na;
p_in(:,3)=p_in(:,3)-.002; 
p_in=p_in(:,[1 3 2]);
 tr_in=tr_in(:,[1,3,2]);
p_in(:,3)=p_in(:,3); 
nin=numel(tr_in)/3;

%combine axon and electrode arrays
tri=cat(1,tr_in,tr_dbs+max(tr_in(:)))';
p=cat(1,p_in,p_dbs)';
elecid2=cat(1,zeros(size(tr_in(:,1))),elecid_dbs);

%plot the geometry of this scenario
trisurf(tri',p(1,:)',p(2,:)',p(3,:)',elecid2)

%Generated the DBS electrode right-hand side
RHS=DBS_RHS(elecid2,tri,p);



%% Define material parameters
rl=1; %resistivity of the cytoplasm
sig=[1/rl 2]; %conductivity array
sigav=cat(1,(sig(1)+sig(2))/2*ones([nin 1]),(0+sig(2))/2*ones([ndbs 1])); % average conductivity accross boundary array
sigdiff=cat(1,(sig(2)-sig(1))*ones([nin 1]),(sig(2)-0)*ones([ndbs 1])); % difference in conductivity accross boundary array
rhotoIm=cat(1,(1/sig(2)-1/sig(1))^-1*ones([numel(p_in)/3 1]),zeros([numel(p_dbs)/3 1]));   %conversion factor for rho to Im


%% Solve the matrix equation at every time step
% Solve bidomain BEM equations and find the threshold
tic
[ Vm Im rho t_mem p_mem tria pa BEM]=bidomainBEM(RHS,wave,1,tarr,'CN2',p,tri,sigav,sigdiff,rhotoIm,1:(numel(p_in)/3),sig);
toc

%construct helper functions
getRow = @(data, rowNum) reshape(data(:,1000:N),[],1); % Helper function
thfunc=@(amp) max(getRow(bidomainBEM(RHS,wave,amp,tarr,'CN2',p,tri,sigav,sigdiff,rhotoIm,1:numel(p_in)/3,sig,BEM)))>0;

%find threshold
BEMampval=findthreshold(@(amp) thfunc(amp),0.020,0.010,.01)

%solve at threshold
tic
[ Vm Im rho t_mem p_mem tria pa BEM]=bidomainBEM(RHS,wave,BEMampval,tarr,'CN2',p,tri,sigav,sigdiff,rhotoIm,1:numel(p_in)/3,sig,BEM);
toc

%% Spatial Plots at 1/5th and 4/5th the axon length
tarr = 1000*tarr;
Vm = 1000*Vm;
test_figure = figure;
xlabel("time (ms)",'FontSize', 10)
plot(tarr,Vm(id1,:),'LineStyle','-','Color',"#ADD8E6",'LineWidth', 3) %light blue
hold on
plot(tarr,Vm(id4,:),'LineStyle','-','Color',"#000080",'LineWidth', 3)% blue
ylabel("Voltage (mV)",'FontSize', 10)
ylim([-80 50])
hold off
xlabel("time (ms)",'FontSize', 10)
test_figure.Units='inches';
test_figure.Position=[-8.0104 4.9792 2 2];
lgd2=legend('BEM-1/5','BEM-4/5');
lgd2.NumColumns = 2;