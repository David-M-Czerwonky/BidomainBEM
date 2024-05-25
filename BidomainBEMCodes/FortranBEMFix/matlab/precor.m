function [S0,D0,Ach,Alin,rfmm,nhatfmm]=precor(t2p,p,col,row,qord)
%this function computes matrix entries associated with all integral operators
if numel(t2p(:,1))~=3
t2p=t2p';
p=p';
end
nt=numel(t2p)/3;
np=numel(p)/3;
ncols=numel(col);
v1=p(:,t2p(1,:))-p(:,t2p(3,:));
v2=p(:,t2p(2,:))-p(:,t2p(3,:));
normal=cross(v1,v2,1);
area=sqrt(normal(1,:).^2+normal(2,:).^2+normal(3,:).^2);
normal(1,:)=normal(1,:)./area;
normal(2,:)=normal(2,:)./area;
normal(3,:)=normal(3,:)./area;
area=area/2;
[qwt,qpt,nquad]=integrationrules2d(qord);
qwts0=zeros([nquad nt]);
qwts1=zeros([nquad nt 3]);
rfmm=zeros([3 nquad nt]);
nhatfmm=zeros([3 nquad nt]);
col0=zeros([nquad nt]);
row0=zeros([nquad nt]);
col1=zeros([nquad nt 3]);
row1=zeros([nquad nt 3]);
rowint=zeros([ncols nquad nquad]);
colint=zeros([ncols nquad nquad]);
for i=1:nquad
rfmm(1,i,:)=reshape(...
            qpt(i,1)*p(1,t2p(1,:))+...
            qpt(i,2)*p(1,t2p(2,:))+...
            qpt(i,3)*p(1,t2p(3,:)),[1 1 nt]);
rfmm(2,i,:)=reshape(...
            qpt(i,1)*p(2,t2p(1,:))+...
            qpt(i,2)*p(2,t2p(2,:))+...
            qpt(i,3)*p(2,t2p(3,:)),[1 1 nt]);
rfmm(3,i,:)=reshape(...
            qpt(i,1)*p(3,t2p(1,:))+...
            qpt(i,2)*p(3,t2p(2,:))+...
            qpt(i,3)*p(3,t2p(3,:)),[1 1 nt]);
nhatfmm(:,i,:)=reshape(normal,[3,1,nt]);
col0(i,:)=i+(0:nt-1)*nquad;
row0(i,:)=(1:nt);
qwts0(i,:)=area*qwt(i);
col1(i,:,1)=i+(0:nt-1)*nquad;
col1(i,:,2)=i+(0:nt-1)*nquad;
col1(i,:,3)=i+(0:nt-1)*nquad;
row1(i,:,1)=1+3*(0:nt-1);
row1(i,:,2)=2+3*(0:nt-1);
row1(i,:,3)=3+3*(0:nt-1);
qwts1(i,:,1)=qwts0(i,:)*qpt(i,1);
qwts1(i,:,2)=qwts0(i,:)*qpt(i,2);
qwts1(i,:,3)=qwts0(i,:)*qpt(i,3);
for j=1:nquad

colin(:,i,j)=i+(col(:)-1)*nquad;
rowin(:,i,j)=j+(row(:)-1)*nquad;
end
end
Ach=sparse(col0(:),row0(:),qwts0(:),nquad*nt,nt);
Alin=sparse(col1(:),row1(:),qwts1(:),nquad*nt,3*nt);
clear col0 row0 qwts0 col1 row1 qwts1;
nfmm=nt*nquad;
nquad2=nquad*nquad;
mex_id_ = 'precorr(i double[xx], i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], o double[xx], o double[xx])';
[S0, D0] = BEM(mex_id_, rfmm, nt, nquad, row, col, ncols, normal, 3, nfmm, 1, 1, ncols, ncols, 1, 3, nt, ncols, nquad2, ncols, nquad2);
S0=sparse(colin(:),rowin(:),S0(:)/(4*pi),nquad*nt,nquad*nt);
D0=sparse(colin(:),rowin(:),D0(:)/(4*pi),nquad*nt,nquad*nt);
clear colin,rowin;
rfmm=reshape(rfmm,[3,nquad*nt]);
nhatfmm=reshape(nhatfmm,[3,nquad*nt]);
end


