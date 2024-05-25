function rhsF1=DBSRHS0(elecid,tri,p)
eps0=8.8541878128*10^-12;
nt=numel(tri)/3;
np=numel(p)/3;

rhsF1=zeros([nt 1]);
tri=tri';p=p';
nhat=cross(p(tri(:,2),:)-p(tri(:,1),:),...
           p(tri(:,3),:)-p(tri(:,1),:));

ap=0;an=0;
for i=1:nt
if sign(elecid(i))==1
    ap=ap+norm(nhat(i,:))/2;
elseif sign(elecid(i))==-1
    an=an+norm(nhat(i,:))/2;
end
end

for i=1:nt
if sign(elecid(i))==1
    rhsF1(i)=-eps0*elecid(i)/ap;
elseif sign(elecid(i))==-1
    rhsF1(i)=-eps0*elecid(i)/an;
end
end