function [p,tri,id,V,F] = GenMeshAxon2Box(test)
%Generates a cylindrical axon in a heterogeneous extracellular bath

%axon properties
Ncompart=100;
radius=10^-6;
Naround=8;
length=0.004;

%Generate cylindrical axon
[V,F] = cylinder_mesh(radius,Naround,'Height',length,"Stacks",Ncompart,'Caps','on');
V=V(:,[1 3 2]);
F=F(:,[1 3 2]);

if test==1
figure
trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',1,'FaceColor','b')
% axis equal
% light
% lighting gouraud
end

%generate box
vx=[min(V(:,1)) max(V(:,1)) mean(V(:,1))];

vy=[min(V(:,2)) max(V(:,2)) mean(V(:,2))];

vz=[min(V(:,3)) max(V(:,3)) mean(V(:,3))];

[V2,F2] = cube(5);

V2(:,1)=V2(:,1)*(vx(2)-vx(1))*8000;
V2(:,2)=V2(:,2)*(vy(2)-vy(1))*2;
V2(:,3)=V2(:,3)*(vz(2)-vz(1))*8000;
V3 = 2*V2;

V2(:,1)=V2(:,1)+vx(3)-4000*(vx(2)-vx(1));
V2(:,2)=V2(:,2)+vy(3)-1*(vy(2)-vy(1));
V2(:,3)=V2(:,3)+vz(3)-4000*(vz(2)-vz(1));


F3 = F2;
V3(:,1)=V3(:,1)+vx(3)-8000*(vx(2)-vx(1));
V3(:,2)=V3(:,2)+vy(3)-2*(vy(2)-vy(1));
V3(:,3)=V3(:,3)+vz(3)-8000*(vz(2)-vz(1));

if test==1
hold on
trisurf(F2,V2(:,1),V2(:,2),V2(:,3),'edgealpha',0,'facealpha',0.3,'FaceColor','#FF0000')

trisurf(F3,V3(:,1),V3(:,2),V3(:,3),'edgealpha',0,'facealpha',0.3,'FaceColor','#008000')

xlabel('x')
ylabel('y')
zlabel('z')
end

%combine meshes
[p,tri] = cat_meshes(V,F,V2,F2,V3,F3);
ne = numel(tri)/3;
nf1 = numel(F)/3;
nf2 = numel(F2)/3;
nf3 = numel(F3)/3;
%id each boundary
id = ones([ne 1]);
id(1:(nf1)) = 1;
id((nf1+1):(nf1+nf2)) = 2;
id((nf1+nf2+1):(nf1+nf2+nf3)) = 3;
end