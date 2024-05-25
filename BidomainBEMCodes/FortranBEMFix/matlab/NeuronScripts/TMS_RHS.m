function  rhs= TMS_RHS(t2p,p,nhat,sigdif,rs,js,area)
np = numel(p)/3;
nt = numel(t2p)/3;
[qwt,qpt,nquad]=integrationrules2d(1);
pquad=[];
for j=1:nquad
    pquad=cat(1,pquad,qpt(j,1)*p(t2p(:,1),:)+qpt(j,2)*p(t2p(:,2),:)+qpt(j,3)*p(t2p(:,3),:));
end

Eprimary=Eprimfmm(10^-7,rs',js',pquad')';%compute Efield samples at nodes

rhs=zeros(np,1);
for i=1:nt
    for j=1:nquad
        con=sigdif(i)*sum(Eprimary(i+(j-1)*nt,:).*nhat(i,:))*area(i)*qwt(j);
        rhs(t2p(i,1))=rhs(t2p(i,1))+con*qpt(j,1);
        rhs(t2p(i,2))=rhs(t2p(i,2))+con*qpt(j,2);
        rhs(t2p(i,3))=rhs(t2p(i,3))+con*qpt(j,3);
    end
end
end

