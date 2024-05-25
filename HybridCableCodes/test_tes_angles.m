%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Boundary Element Method of Bidomain Modeling for Predicting Cellular
%                 Responses to Electromagnetic Fields
%
%             Code Author: David M. Czerwonky
%
% Description: This script runs the hybrid cable approach to predict the
% behavior of HH axon embeded in a rectangular extracellular bath that is
% exposed to TES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('./Meshing/')
addpath('./Meshing/matrix')
addpath('./Meshing/mesh/')

%load axon properties
a=10^-6;   %axon radius
len=.004;  %axon length
rl=1;      %axon longitudinal resistivity
M=100;     %number of compartments

%time scale
dt=1*10^-6;   %time step
dur=20*10^-3;   % duration

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

%find the actication threshold (thval) and the transmembrane voltage (Vm)
%over a full 90 deg axon rotation
thval=tes_angles(a,len,rl,Iinject,phead',t2phead',elecid,epseff,dt,dur,M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                             End of Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%