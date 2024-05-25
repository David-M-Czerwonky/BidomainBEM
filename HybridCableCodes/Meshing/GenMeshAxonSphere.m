function [p,tri,id,V,F] = GenMeshAxonSphere(test,res)
%Generates the a cylindrical axon in a spherical bath

%axon properties
Ncompart=100;
radius=10^-6;
Naround= 8;
length=.004;

%Generate Axon
[V,F] = cylinder_mesh(radius,Naround,'Height',length,"Stacks",Ncompart,'Caps','on');
V=V(:,[1 3 2]) + [0 -0.002 0.07];
% V = rotate(rotate(V,pi/2,"y")',pi/2,"x")';
F=F(:,[1 3 2]);

if test==1
figure
trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',1)
hold on
% axis equal
% light
% lighting gouraud
end

%Generate Extracellular Space as a Sphere
[V2,F2] = icosphere(res);
V2 = 0.085*V2; %make it an 85mm sphere which is commonly used in TMS studies
if test==1
trisurf(F2,V2(:,1),V2(:,2),V2(:,3),'edgealpha',0,'facealpha',0.3,'FaceColor','#FF0000')
xlabel('x')
ylabel('y')
zlabel('z')
end

%combine meshes
[p,tri] = cat_meshes(V,F,V2,F2);
ne = numel(tri)/3;
nf1 = numel(F)/3;
nf2 = numel(F2)/3;
%id each boundary
id = ones([ne 1]);
id(1:(nf1)) = 1;
id((nf1+1):(nf1+nf2)) = 2;
end