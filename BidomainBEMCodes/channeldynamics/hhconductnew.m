function nmh=hhconductnew(Vm,nmhp,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Boundary Element Method of Bidomain Modeling for Predicting Cellular
%                 Responses to Electromagnetic Fields
%
%             Code Author: David M. Czerwonky
%
% Description: This function solves the auxilliary gating equations using
% the exponential Euler rule on pg 193 of Theortical Neuroscience by Dayan
% and Abbott 2001
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt=1000*dt;
Vm=1000*Vm;

%calculate the ion channel rate functions
alphaN= 0.01*(Vm+55)./(1-exp(-0.1*(Vm+55)));
alphaM= 0.1*(Vm+40)./(1-exp(-0.1*(Vm+40)));
alphaH= 0.07*exp(-0.05*(Vm+65));
betaN=  0.125*exp(-0.0125*(Vm+65));
betaM=  4*exp(-0.0556*(Vm+65));
betaH= 1./(1+exp(-.1*(Vm+35)));

%interpolate the NaN values so the function is continuous
alphaN(isnan(alphaN))=.1;
alphaM(isnan(alphaM))=1;
%if any rate function tends to infinity, we cap the value at 10^90 to avoid
% NaNs and infs in the time and steady state limit
alphaN(isinf(alphaN))=10^90;
alphaM(isinf(alphaM))=10^90;
alphaH(isinf(alphaH))=10^90;
betaN(isinf(betaN))=10^90;
betaM(isinf(betaM))=10^90;
betaH(isinf(betaH))=10^90;

%calculate the gating time "constants"
tauN= 1./(alphaN+betaN);
tauM= 1./(alphaM+betaM);
tauH= 1./(alphaH+betaH);

%calculate the steady state limit of the gating variable
Ninfti= alphaN.*tauN;
Minfti= alphaM.*tauM;
Hinfti= alphaH.*tauH;

%Solve the equations using the exponential Euler update rule
nmh=cat(2,...
    Ninfti+(nmhp(:,1)-Ninfti).*exp(-dt./tauN),...
    Minfti+(nmhp(:,2)-Minfti).*exp(-dt./tauM),...
    Hinfti+(nmhp(:,3)-Hinfti).*exp(-dt./tauH));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                             End of Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
