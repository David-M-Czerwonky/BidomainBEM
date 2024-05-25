%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Boundary Element Method of Bidomain Modeling for Predicting Cellular
%                 Responses to Electromagnetic Fields
%
%             Code Author: David M. Czerwonky
%
% Description: This script runs a simulatio for TES for an axon in an
% extracellular bath as shown in figure 5 of the publication. Less data
% points are included here and we use guess values from a similar
% simulation to reduce the time it takes to find the new thresholds.
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
amp = 11.6*10^-3;
%waveform information
close all
Iinject=@(t) (t/(5*10^-5).*(t<5*10^-5)+(t>=5*10^-5)).*(t<5*10^-3);
wave=zeros(size(tarr));
wave(1:numel(tarr))=Iinject(tarr-0.0002).*(tarr>0.0002);

%load some guess values based on hybrid cable results for a similar
%waveform (thval)
load('TESAngleGuess.mat')

BEMthval = zeros(19,1);
angle = 0:5:90;

for j = 1:19
    %% Initalize Geometry
    mesh_print = 0;
    res = 10;
    [p,tri,id,V,F] = GenMeshAxonAngleBox(mesh_print,res,angle(j)/180*pi);
    
    %% Define material parameters
    rl=1; %resistivity of the cytoplasm
    sig=[1/rl 2];
    sigav=zeros(numel(id),1);
    sigdiff=zeros(numel(id),1);
    rhotoIm=zeros(numel(p)/3,1);
    
    sigav(id==1) = (sig(1)+sig(2))/2;
    sigav(id==2) = (sig(2)+0)/2;
    
    sigdiff(id==1) = (sig(2)-sig(1));
    sigdiff(id==2) = (0-sig(2));
    
    Nneu = numel(V)/3; %number of nodes on the neuron membrane
    rhotoIm(1:Nneu) = 1./(1./sig(2) -1./sig(1));
    rhotoIm((Nneu+1):end) = 0;
    
    %Generate TES electrode forcing function
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
    RHS=TES_RHS(elecid2,tri',p');
    
    %solve bidomain BEM according to the flow chart in publication to
    %initalize the BEM matrix data structures
    [ Vmem2 Im2 rho2 t_mem2 p_mem1 tria pa BEM2]=bidomainBEM(RHS,wave,-1,tarr,'CN2',p',tri',sigav,sigdiff,rhotoIm,1:(numel(V)/3),sig);
    
    %created helper functions to to find the threshold
    getRow = @(data, rowNum) reshape(data(:,1000:N),[],1); % Helper function
    thfunc=@(amp) max(getRow(bidomainBEM(RHS,wave,-amp,tarr,'CN2',p',tri',sigav,sigdiff,rhotoIm,1:numel(V)/3,sig,BEM2)))>0;
    
    %use the bisection method to find the threshold predicted by bidomain
    %BEM for each angle.
    BEMthval(j)=findthreshold(@(amp) thfunc(amp),100*thval(j),0.01*thval(j),.01)
end