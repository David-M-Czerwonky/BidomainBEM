function [p,tri]=cylmesh(Na,Nl,len,wid,test)
[x,y]=ndgrid(0:2*pi/Na:2*pi-2*pi/Na,0:1/(Nl-1):1);
tri=delaunay(x(:),y(:));

ival=0:Nl-2;

tri=cat(1,tri,[Na*(ival(:)+1),Na*(ival(:))+1,Na*(ival(:)+2)],[Na*(ival(:))+1,Na*(ival(:)+1)+1,Na*(ival(:)+2)]);
tri=cat(1,tri,[([2:Na,1])',(1:Na)',(Na*Nl+1)*ones([Na 1])]);

tri=cat(1,tri,[(1:Na)'+Na*(Nl-1),([2:Na,1])'+Na*(Nl-1),(Na*Nl+2)*ones([Na 1])]);
clear p
p(:,1)=wid*cos(x(:));
p(:,2)=wid*sin(x(:));
p(:,3)=len*y(:);
p(end+1,:)=[0 0 0];
p(end+1,:)=[0 0 len];

if test==1
hold on
trisurf(tri,p(:,1),p(:,2),p(:,3))
p=p';tri=tri';
v1=p(:,tri(1,:))-p(:,tri(3,:));
v2=p(:,tri(2,:))-p(:,tri(3,:));
normal=cross(v1,v2,1)';
pcen=(p(:,tri(1,:))+p(:,tri(2,:))+p(:,tri(3,:)))'/3;
quiver3(pcen(:,1),pcen(:,2),pcen(:,3),normal(:,1),normal(:,2),normal(:,3))
p=p';tri=tri';
end