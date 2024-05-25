function rhsF1=DBS_RHS(elecid,tri,p)
nt=numel(tri)/3;
np=numel(p)/3;

rhsF1=zeros([np 1]);
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
        rhsF1(tri(i,:))=rhsF1(tri(i,:))-norm(nhat(i,:))/6/ap*elecid(i);
    elseif sign(elecid(i))==-1
        rhsF1(tri(i,:))=rhsF1(tri(i,:))-norm(nhat(i,:))/6/an*elecid(i);
    end
end