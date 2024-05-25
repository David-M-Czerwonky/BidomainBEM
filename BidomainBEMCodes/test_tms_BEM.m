%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Boundary Element Method of Bidomain Modeling for Predicting Cellular
%                 Responses to Electromagnetic Fields
%
%             Code Author: David M. Czerwonky
%
% Description: This script simulates a cylindrical axon immersed in a
% spherical extracellular bath and stimlated by a coil 5cm above the bath
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
dt=10^-6;     %time step
dur=20*10^-3; % duration
tarr=0:dt:dur;
N=numel(tarr);

%% Initalize the TMS coil model
%load coil source information
load coilinformation.mat;
rv(:,3)=rv(:,3)+.09; 

%initalize the device waveform
waveform=@(t) VL(3000,t-.0002)/L.*(t-.0002>0);
wave = waveform(tarr);

%% Initalize Geometry
mesh_print = 0; %This will print a figure of the geometry if 1, otherwise it will not print
res = 3;        %This decides the resolution of the spherical bath
[p,tri,id,V,F] = GenMeshAxonSphere(mesh_print,res);

%% Define material parameters
rl=1; %resistivity of the cytoplasm
sig=[1/rl 2]; %conductivity array
sigav=zeros(numel(id),1); %average conductivity across boundary array
sigdiff=zeros(numel(id),1); %conductivity difference across boundary array
rhotoIm=zeros(numel(p)/3,1); %conversion factor for rho to Im

sigav(id==1) = (sig(1)+sig(2))/2;
sigav(id==2) = (sig(2)+0)/2;

sigdiff(id==1) = (sig(2)-sig(1));
sigdiff(id==2) = (0-sig(2));

Nneu = numel(V)/3; %number of nodes on the neuron membrane
rhotoIm(1:Nneu) = 1./(1./sig(2) -1./sig(1));
rhotoIm((Nneu+1):end) = 0;

%% Construct the TMS right-hand side of the matrix equation
p=p'; tri=tri';
%construct the boundary normal vectors and element areas
v1=p(:,tri(1,:))-p(:,tri(3,:));
v2=p(:,tri(2,:))-p(:,tri(3,:));
normal=cross(v1,v2,1);
area=sqrt(normal(1,:).^2+normal(2,:).^2+normal(3,:).^2);
normal(1,:)=normal(1,:)./area;
normal(2,:)=normal(2,:)./area;
normal(3,:)=normal(3,:)./area;
area = 0.5*area;
p=p'; tri=tri'; normal = normal';
%Construct right-hand side for TMS coil
RHS=TMS_RHS(tri,p,normal,sigdiff,rv,Jdc,area);

%% Solve the matrix equation at every time step
%get run time for one complete solution across the entire duration
tic
[ Vm Im rho t_mem p_mem tria pa BEM]=bidomainBEM(RHS,wave,-1,tarr,'ECN',p',tri',sigav,sigdiff,rhotoIm,1:(numel(V)/3),sig);
toc

%construct helper functions
getRow = @(data, rowNum) reshape(data(:,1000:N),[],1); % Helper function
thfunc=@(amp) max(getRow(bidomainBEM(RHS,wave,-amp,tarr,'ECN',p',tri',sigav,sigdiff,rhotoIm,1:numel(V)/3,sig,BEM)))>0;

%find the threshold predicted by bidomain BEM
BEMampval=findthreshold(@(amp) thfunc(amp),7,5,.01)

%solve at threshold
tic
[ Vm Im rho t_mem p_mem tria pa BEM]=bidomainBEM(RHS,wave,-BEMampval,tarr,'ECN',p',tri',sigav,sigdiff,rhotoIm,1:numel(V)/3,sig,BEM);
toc
max(max(Vm))

%% Spatial Plots at 1/5th and 4/5th the axon length
tarr = 1000*tarr;
Vm = 1000*Vm;
test_figure = figure;
id15 = ceil(1/5*810);
id45 = ceil(4/5*810);
xlabel("time (ms)",'FontSize', 10)
plot(tarr,Vm(id15,:),'LineStyle','-','Color',"#ADD8E6",'LineWidth', 3) %light blue
hold on
plot(tarr,Vm(id45,:),'LineStyle','-','Color',"#000080",'LineWidth', 3)% blue
ylabel("Voltage (mV)",'FontSize', 10)
ylim([-80 50])
hold off
xlabel("time (ms)",'FontSize', 10)
test_figure.Units='inches';
test_figure.Position=[-8.0104 4.9792 2 2];
lgd2=legend('BEM-1/5','BEM-4/5');
lgd2.NumColumns = 2;