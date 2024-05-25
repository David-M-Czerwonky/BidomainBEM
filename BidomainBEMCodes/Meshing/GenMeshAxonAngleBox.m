function [p,tri,id,V,F] = GenMeshAxonAngleBox(test,res,th)
%Generates a cylindrical axon in an extracellular bath at a certain
%orientation

%axon properties
Ncompart=100;
radius=10^-6;
Naround=8;
length=0.004;

%Generate cylindrical axon
[V,F] = cylinder_mesh(radius,Naround,'Height',length,"Stacks",Ncompart,'Caps','on');
V=V(:,[1 3 2]) -[0 2*10^-3 0];
F=F(:,[1 3 2]);
ptemp = V;
V(:,1)=ptemp(:,1);
V(:,2)=ptemp(:,2)*cos(th)-ptemp(:,3)*sin(th);
V(:,3)=ptemp(:,2)*sin(th)+ptemp(:,3)*cos(th);

if test==1
figure
trisurf(F,V(:,1),V(:,2),V(:,3),'edgealpha',1,'FaceColor','b')
end

%generate box
[V2,F2,~] = cube(res);
V2(:,1)=V2(:,1)-1/2;
V2(:,2)=V2(:,2)-1/2;
V2(:,3)=V2(:,3)-1/2;
V2(:,1)=.016*V2(:,1);
V2(:,2)=.008*V2(:,2);
V2(:,3)=.016*V2(:,3);

%combine meshes
[p,tri] = cat_meshes(V,F,V2,F2);
ne = numel(tri)/3;
nf2 = numel(F2)/3;
id = 2*ones([ne 1]);
id(1:(ne-nf2)) = 1;

if test==1
hold on
trisurf(F2,V2(:,1),V2(:,2),V2(:,3),'edgealpha',0,'facealpha',0.3,'FaceColor',[150 75 0]/256)
xlabel('x')
ylabel('y')
zlabel('z')
view(90,0)
end

end