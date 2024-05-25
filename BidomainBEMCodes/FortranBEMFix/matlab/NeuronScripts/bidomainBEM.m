function [ Vmem Imstore rhostore t_mem p_mem tri p BEM2 rhoRHSeul rhoRHSAmat]=bidomainBEM(RHS,wave,amp,tarr,type,p,tri,sigav,sigdiff,rhotoIm,p_mem,sig,BEM2)

%% Variable Initalizations
Ntime=numel(tarr);
len=sqrt(sum((p(:,tri(1,:))-p(:,tri(2,:))).^2,1));
np=numel(p(:))/3;
nbdr=numel(p_mem);
conv=sparse(1:np,1:np,rhotoIm);
t_mem=tri(:,abs(sigav-(mean(sig(1:2))))<10^-8);
[p_mem,~,t_mem]=unique(t_mem(:));
t_mem=reshape(t_mem,[3 numel(t_mem)/3]);
nt=numel(tri)/3;
np=numel(p)/3;
fmmquad=6;
nearlen=10^16;

%%% HH membrane parameters
cm=10^-2;
gl=3;gk=360;gna=1200;
El=-0.0544;Ek=-0.077;Ena=0.05;

%%%initial membrane state conditions
nmhp=zeros([nbdr 3]);
nmhp(:,1)=3.176862580749554e-01;
nmhp(:,2)=5.294913493474710e-02; 
nmhp(:,3)=5.963697918631813e-01;

%%%% Compute BEM2 matrices
dt=tarr(2)-tarr(1);
c=dt./(2*cm);
if ~exist('BEM2','var')
    % BEM2 parameter does not exist, so it must be calculated
    BEM2=adjneursetup2(tri,p,sigdiff,sigav,nearlen,fmmquad,c,conv,p_mem);
    [u,s,v]=svd(BEM2.Amat);
    s=diag(diag(1./s));s(end)=0;
    BEM2.Amat=v*s*u';
    [u,s,v]=svd(BEM2.AnearD);
    s=diag(diag(1./s));s(end)=0;
    BEM2.AnearD=v*s*u';
end
b=amp*RHS;
rhoRHSAmat=BEM2.Amat*(BEM2.preconAmat*b);
rhoRHSeul=BEM2.AnearD*(BEM2.preconeul*b);

%Initialize all boundary quantities that we want to solve for
Im=zeros([nbdr 1]);
Vmem=zeros([nbdr Ntime]);
Imstore=zeros([nbdr Ntime]);
Vmem(:,1)=-.065;
Vm=Vmem(:,1);
rho=zeros([np 1]);
rhostore=zeros([np Ntime]);

%% Time Integration
for t=1:Ntime-1
    nmh=nmhp;
    giter=(gl+gk*nmh(:,1).^4+gna*nmh(:,2).^3.*nmh(:,3));
    Ein=(gl*El+gk*nmh(:,1).^4*Ek+gna*nmh(:,2).^3.*nmh(:,3)*Ena);
    dt=tarr(t+1)-tarr(t);
    c=dt/(2*cm);
    
    if type=='ECN'
        Iion=giter.*Vm-Ein;
        g=Vm-2*c*(Iion);
        [Im,rho]=findImeul2(BEM2,BEM2.AnearD,wave(t),Vm,conv,p_mem,rho,rhoRHSeul*wave(t));
        Im=Im(p_mem);
        Imstore(:,t)=Im;
        g=g+c*Im;
        [Im,rho]=findImnicolsson2(BEM2,BEM2.Amat,wave(t+1)*b,g,conv,p_mem,rho,rhoRHSAmat*wave(t+1));
        Im=Im(p_mem);
        Vm=c*Im+g;
        Vmem(:,t+1)=Vm;
    elseif type=='CN2'
        Iion=giter.*Vm-Ein;
        g=Vm+dt/cm*(.5*Im-Iion);
        [Im,rho]=findImnicolsson2(BEM2,BEM2.Amat,wave(t+1)*b,g,conv,p_mem,rho,rhoRHSAmat*wave(t+1));
        Im=Im(p_mem);
        Vm=0.5*dt/cm*Im+g;
        Vmem(:,t+1)=Vm;
        Imstore(:,t+1)=Im;
        rhostore(:,t+1)=rho(:);
    elseif type=='EUL'
        dt=tarr(t+1)-tarr(t);
        [Im,rho]=findImeul2(BEM2,BEM2.AnearD,wave(t),Vm,conv,p_mem,rho,rhoRHSeul*wave(t));
        Im=Im(p_mem);
        Imstore(:,t)=Im;
        Iion=giter.*Vm-Ein;
        %euler
        Vm=Vm+dt/cm.*(Im-Iion);
        Vmem(:,t+1)=Vm;
    end
    % Update Gating Variables
    nmhp=hhconductnew(Vm,nmhp,dt);

end
% BEM_activity = max(max(Vmem))
end