function lap=grammsurf(t2p,p,sigav)
np=numel(p)/3;
nt=numel(t2p)/3;
col=ones([nt*9 1]);
row=ones([nt*9 1]);
val=zeros([nt*9 1]);
for i=1:nt
col(9*(i-1)+1:9*i)=t2p(i,[1 1 1  2 2 2  3 3 3 ]);
row(9*(i-1)+1:9*i)=t2p(i,[1 2 3  1 2 3  1 2 3 ]);
Area=0.04166666666666666*norm(cross(...
p(t2p(i,2),:)-p(t2p(i,1),:),p(t2p(i,3),:)-p(t2p(i,1),:)));

val(9*(i-1)+1:9*i)=sigav(i)*[2,1,1,1,2,1,1,1,2]*Area;
end
lap=sparse(col,row,val,np,np);