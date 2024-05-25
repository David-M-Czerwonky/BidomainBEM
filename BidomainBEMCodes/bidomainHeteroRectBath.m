%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Boundary Element Method of Bidomain Modeling for Predicting Cellular
%                 Responses to Electromagnetic Fields
%
%             Code Author: David M. Czerwonky
%
% Description: This script predicts E-fields without including a membrane
% boundary. We find the E-field predicted by bidomain BEM and the standard
% Adjoint BEM. 
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
dur1=10*10^-2;
dur=20*10^-3;
dt=4*10^-6;
tarr=0:dt:dur;
N=numel(tarr);

%% Initalize the TES model
%waveform information
close all
Iinject=@(t) (t/(5*10^-5).*(t<5*10^-5)+(t>=5*10^-5)).*(t<5*10^-3);
wave=zeros(size(tarr));
wave(1:numel(tarr))=Iinject(tarr-0.0002).*(tarr>0.0002);


%% Initalize Geometry
load('2LayerRectBath.mat','p','tri','id','V','F')
% trisurf(tri,p(:,1),p(:,2),p(:,3),'edgealpha',0.2,'facealpha',0.2,'facecolor','#FF0000')


%% Define material parameters
rl=1; %resistivity of the cytoplasm
sig=[1/rl 0.276 1.654]; %conductivity array
sigav=zeros(numel(id),1); %average conductivity across boundary array
sigdiff=zeros(numel(id),1); %conductivity difference across boundary array
rhotoIm=zeros(numel(p)/3,1); %conversion factor for rho to Im

sigav(id==1) = (sig(1)+sig(2))/2;
sigav(id==2) = (sig(2)+sig(3))/2;
sigav(id==3) = (sig(3)+0)/2;

sigdiff(id==1) = (sig(2)-sig(1));
sigdiff(id==2) = (sig(3)-sig(2));
sigdiff(id==3) = (0-sig(3));

Nneu = numel(V)/3; %number of nodes on the neuron membrane
rhotoIm(1:Nneu) = 1./(1./sig(2) -1./sig(1));
rhotoIm((Nneu+1):end) = 0;

%% Construct the TES right-hand side of the matrix equation
%define electrode locations
elecid2 = zeros([numel(tri)/3 1]);
pcen=(p(tri(:,1),:)+p(tri(:,2),:)+p(tri(:,3),:))/3;
elecdir=2;
for i=1:numel(pcen)/3
    if(pcen(i,elecdir)<=min(pcen(:,elecdir))+.000001)
        elecid2(i) = 1;
    elseif(pcen(i,elecdir)>=max(pcen(:,elecdir))-.000001)
        elecid2(i) = -1;
    end
end
%construct TES right-hand side
RHS=TES_RHS(elecid2,tri',p');


%% Solve the matrix equation at every time step
%get run time for one complete solution across the entire duration
tic
[ Vm Im rho t_mem p_mem tria pa BEM]=bidomainBEM(RHS,wave,-1,tarr,'CN2',p',tri',sigav,sigdiff,rhotoIm,1:(numel(V)/3),sig);
toc

%construct helper functions
getRow = @(data, rowNum) reshape(data(:,10:N),[],1); % Helper function
thfunc=@(amp) max(getRow(bidomainBEM(RHS,wave,-amp,tarr,'CN2',p',tri',sigav,sigdiff,rhotoIm,1:numel(V)/3,sig,BEM)))>0;

%find the threshold predicted by bidomain BEM
BEMampval=findthreshold(@(amp) thfunc(amp),0.02,0.015,.01)

%solve at threshold
tic
[ Vm Im rho t_mem p_mem tria pa BEM]=bidomainBEM(RHS,wave,BEMampval,tarr,'CN2',p',tri',sigav,sigdiff,rhotoIm,1:numel(V)/3,sig,BEM);
toc
max(max(Vm))

%Spatial Plots
tarr = 1000*tarr;
Vm = 1000*Vm;
test_figure = figure;
id15 = ceil(1/5*810);
id45 = ceil(4/5*810);
xlabel("time (ms)",'FontSize', 10)
plot(tarr,Vm(id15,:),'LineStyle','-','Color',"#ADD8E6",'LineWidth', 3) %light blue
hold on
plot(tarr,Vm(id45,:),'LineStyle','-','Color',"#000080",'LineWidth', 3)% blue
% xlim([0,16])
ylabel("Voltage (mV)",'FontSize', 10)
ylim([-80 50])
hold off
xlabel("time (ms)",'FontSize', 10)
test_figure.Units='inches';
test_figure.Position=[-8.0104 4.9792 2 2];
lgd2=legend('BEM-1/5','BEM-4/5');
lgd2.NumColumns = 2;