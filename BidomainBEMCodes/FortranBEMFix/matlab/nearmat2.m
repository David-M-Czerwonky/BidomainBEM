function [S0,Dstar0,S1,D1,col,row,area]=nearmat2(t2p,p,dnear)
%this function computes matrix entries associated with all integral operators
if numel(t2p(:,1))~=3
t2p=t2p';
p=p';
end
nt=numel(t2p)/3;
np=numel(p)/3;
v1=p(:,t2p(1,:))-p(:,t2p(3,:));
v2=p(:,t2p(2,:))-p(:,t2p(3,:));
normal=cross(v1,v2,1);
area=sqrt(normal(1,:).^2+normal(2,:).^2+normal(3,:).^2);
normal(1,:)=normal(1,:)./area;
normal(2,:)=normal(2,:)./area;
normal(3,:)=normal(3,:)./area;
area=area/2;
del=dnear*sqrt(mean(area(:))*2);
ncol=nt^2;
[row col]=ndgrid(1:nt,1:nt);
row=row(:);col=col(:);
row=row(col~=0);
col=col(col~=0);
ncols=numel(col);
lincol=zeros(ncols,9);
linrow=zeros(ncols,9);
for i=1:3
for j=1:3
lincol(:,i+(j-1)*3)=3*(col(:)-1)+i;
linrow(:,i+(j-1)*3)=3*(row(:)-1)+j;
end
end
Dstar0=zeros([ncols 1]);
S0=zeros([ncols 1]);
S1=zeros([ncols 9]);
D1=zeros([ncols 9]);
mex_id_ = 'nearmat(i double[xx], i int[x], i int[xx], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[x], o double[x], o double[x], o double[xx], o double[xx])';
[Dstar0, S0, D1, S1] = BEM(mex_id_, p, np, t2p, nt, row, col, ncols, normal, area, 3, np, 1, 3, nt, 1, ncols, ncols, 1, 3, nt, nt, ncols, ncols, ncols, 9, ncols, 9);
Dstar0=full(sparse(col,row,Dstar0/(4*pi)));
S0=full(sparse(col,row,S0/(4*pi)));
D1=full(sparse(lincol(:),linrow(:),D1(:)/(4*pi)));
S1=full(sparse(lincol(:),linrow(:),S1(:)/(4*pi)));
end


