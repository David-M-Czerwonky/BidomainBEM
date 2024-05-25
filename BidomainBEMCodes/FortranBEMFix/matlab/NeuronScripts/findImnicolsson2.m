function [Im,rhot]=findImnicolsson2(BEM,Amat,RHS,g,Conversion,pmem,rhop,rhop2)

V=BEM.preconAmat*(-BEM.AnearS0(:,pmem)*g);
rho=Amat*V; %rho from the hyper singular forcing term
%rhop2 is from the device forcing terms
%By superposition add rho + rhop2 to get the total rho
rhot = rho+rhop2;
Im=Conversion*(rhot);
end
