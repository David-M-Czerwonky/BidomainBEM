function [thval Vm wave tarr]=tes_angles(a,len,rl,Iinject,phead,t2phead,elecid,epseff,dt,dur,M)
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

%calculate the right-hand side from the DBS electrode
rhs=TESRHS0(elecid,t2phead,phead);

%define the locations of the axon compartments
p=zeros([M 3]);
p(:,3)=-len/2:len/(M-1):len/2;
p=p(:,[1 3 2]);
p(:,3)=p(:,3);

%define the helper functions for findthreshold
getRow = @(data, rowNum) data(rowNum,1000:N); % Helper function


%find the actication threshold (thval) and the transmembrane voltage (Vm)
thval = zeros(91,1);
angle = 0:90;
ptemp=zeros([M 3]);
Ve=zeros([M 1]); % initalize extracellular potential
for i = 90
    %rotate the locations of the axon compartments
    th = angle(i)/180*pi;
    ptemp(:,1)=p(:,1);
    ptemp(:,2)=p(:,2)*cos(th)-p(:,3)*sin(th);
    ptemp(:,3)=p(:,2)*sin(th)+p(:,3)*cos(th);
    
    dl=ptemp(2:end,:)-ptemp(1:end-1,:);
    peval=(ptemp(2:end,:)+ptemp(1:end-1,:))/2;
    %calculate the E-field produced by the DBS electrode
    [Efield]=adjointBEM_electrode(t2phead,phead,epseff,rhs,peval');
    %calculate the extracellular potentials
    Ve(1)=0;
    Ve(2:end)=cumsum(sum(dl'.*Efield,1));
    thfunc=@(amp) max(getRow(multicompartsolv(a,len,rl,M,-amp,Ve,waveform,dt,dur),round(M/5)))>0;
    thval(i)=findthreshold(@(amp,Ve) thfunc(amp),100,10^-3,.01);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                             End of Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%