function rhs=TESRHS(elecid,tri,p)
nt=numel(tri)/3;
np=numel(p)/3;

rhs=zeros([np 1]);
nhat=cross(p(tri(:,2),:)-p(tri(:,1),:),...
           p(tri(:,3),:)-p(tri(:,1),:));
%Find total electrode area
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
        rhs(tri(i,:))=rhs(tri(i,:))-norm(nhat(i,:))/6*elecid(i)/ap;
    elseif sign(elecid(i))==-1
        rhs(tri(i,:))=rhs(tri(i,:))-norm(nhat(i,:))/6*elecid(i)/an;
    end
end
end