function [Im,rho]=findImeul2(BEM,Amat,RHS,g,Conversion,pmem,rhop,rho2)
V=BEM.preconeul*(BEM.AnearS0(:,pmem)*(-g));
rho=Amat*V;
Im=Conversion*(rho+rho2);
end
