function [thval Vm wave tarr]=tes_head(a,len,rl,Iinject,phead,t2phead,elecid,epseff,dt,dur,M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Boundary Element Method of Bidomain Modeling for Predicting Cellular
%                 Responses to Electromagnetic Fields
%
%             Code Author: David M. Czerwonky
%
% Description: This script runs the hybrid cable approach for TES of an
% axon in an extracellular bath.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('./BEM_MEX_Fortran_code')
addpath('./channeldynamics/')
addpath('./Cableequation/')

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
rhs=TESRHS0(elecid,t2phead,phead);

%calculate the E-field produced by the DBS electrode
[Efield]=adjointBEM_electrode(t2phead,phead,epseff,rhs,peval');

%calculate the extracellular potentials
Ve=zeros([M 1]);
Ve(1)=0;
Ve(2:end)=cumsum(sum(dl'.*Efield,1));

%define the helper functions for findthreshold
getRow = @(data, rowNum) data(rowNum,1000:N); % Helper function
thfunc=@(amp) max(getRow(multicompartsolv(a,len,rl,M,-amp,Ve,waveform,dt,dur),round(M/5)))>0;

%find the actication threshold (thval) and the transmembrane voltage (Vm)
thval=findthreshold(@(amp) thfunc(amp),0.05,10^-6,.01);
Vm=multicompartsolv(a,len,rl,M,-thval,Ve,waveform,dt,dur);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                             End of Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%