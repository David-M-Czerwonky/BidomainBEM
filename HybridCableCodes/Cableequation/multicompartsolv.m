function Vm=multicompartsolv(a,len,rl,M,Ve,amp,waveform,dt,dur)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Boundary Element Method of Bidomain Modeling for Predicting Cellular
%                 Responses to Electromagnetic Fields
%
%             Code Author: David M. Czerwonky
%
% Description: This function solves the multi-compartment equations of the
% cable equation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initatlize Time Scale
tarr=0:dt:dur;
N=numel(tarr);

%Membrane properties
cm=10^-2;
gl=3;gk=360;gna=1200;
El=-0.0544;Ek=-0.077;Ena=0.05;
Erest=-0.065;

%Cable length-area parameters
dx=len/(M-1);      %compartment length
AreaL=(pi*a^2)/dx; %longitudinal factor
AreaM=(2*pi*a*dx); %accross membrane factor
AreaMbound=(AreaM+pi*a^2); %end fiber factor
areas=AreaM*ones([M 1]);
areas([1 end])=AreaMbound;

Cm=cm*areas;   %calculate membrane capacitance
Rl=rl/(AreaL); %calculate longitudinal resistance

%create the multi-compartment matrix
A=sparse([1:M,2:M,1:M-1],...
       [1:M,1:M-1,2:M], ...
       [(-2/(Rl))*ones(1,M), ...
        1/Rl*ones(1,M-1), ...
        1/Rl*ones(1,M-1)]);
A(1,1)=(-1/(Rl));
A(end,end)=(-1/(Rl));

%Get activation function from extracellular potentials
Ve=Ve*amp;
Amplitude=A*Ve;
A=.5*A;

%Initalize Vm and auxillary gating variables
Vm=zeros([M N]);
nmhp=zeros([M 3]);
Vm(:,1)=Erest;
nmhp(:,1)=3.176862580749554e-01;
nmhp(:,2)=5.294913493474710e-02; 
nmhp(:,3)=5.963697918631813e-01;

%begin time stepping
for t=1:N-1
    dt=tarr(t+1)-tarr(t);
    nmh=nmhp;
    giter=(gl+gk*nmh(:,1).^4+gna*nmh(:,2).^3.*nmh(:,3)).*areas*0.5;
    Ein=(gl*El+gk*nmh(:,1).^4*Ek+gna*nmh(:,2).^3.*nmh(:,3)*Ena).*areas;
    if dt>10^-7
        %crank-nicolson time stepping
        Vm(:,t+1)=(sparse(1:M,1:M,Cm/dt+giter)-A)\...
          ((Cm/dt-giter).*Vm(:,t)...
          +A*Vm(:,t)+Ein...
          -Amplitude*(waveform(tarr(t+1))+waveform(tarr(t)))/2);
    else
        %forward euler time stepping
         Vm(:,t+1)=Vm(:,t)+dt./Cm.*(...
             -2*giter.*Vm(:,t)+2*A*Vm(:,t)+Ein...
             -Amplitude*waveform(tarr(t)));
    end
    nmhp=hhconductnew(Vm(:,t+1),nmhp,dt);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                             End of Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%