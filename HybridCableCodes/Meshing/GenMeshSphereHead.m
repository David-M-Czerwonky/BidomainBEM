function [p,tri,id] = GenMeshAxonBox(test)
%Generates a cylinder axon in a box
addpath ./neusolv/
addpath('./neusolv/FEM_MEX_C_codes/');
addpath('../')

Ncompart=80;
radius=10^-6;
Naround= 32;

%Generate Cylindrical Axon
[V,F] = cylinder_mesh(radius,Naround,'Height',.004,"Stacks",Ncompart,'Caps','on');
V=V(:,[1 3 2]) + [0 -0.002 0.07];
% V = rotate(rotate(V,pi/2,"y")',pi/2,"x")';
F=F(:,[1 3 2]);

if test==1 || test == 2
figure
trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',1)
hold on
% axis equal
% light
% lighting gouraud
end

%Generate Extracellular Space as a Sphere
[V2,F2] = icosphere(0.085,3);
if test==2
trisurf(F2,V2(:,1),V2(:,2),V2(:,3),'edgealpha',0,'facealpha',0.3)
xlabel('x')
ylabel('y')
zlabel('z')
end

%combine meshes
[p,tri] = cat_meshes(V,F,V2,F2);
ne = numel(tri)/3;
nf2 = numel(F2)/3;
id = ones([ne 1]);
id(1:(ne-nf2)) = -1;

disp("AxonBox Mesh Generation Complete");
end