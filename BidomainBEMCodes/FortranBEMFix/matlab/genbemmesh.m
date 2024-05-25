function [t2p,p2,sigdiff,sigav,rhotoIm]=genbemmesh(te2p,p,condu,regid)
%extracts tissue boundary faces from the tetrahedron mesh
%te2p is the 4 by number of tetrahedrons element array
%p is the 3 by number of points node location array
%condu is the conductivity of reach region
%regid is the number of tetrahedron by one region id array
%t2p is 3 by number of boundary triangles
%p2 is the boundary triangle point array
%sigdiff is the (conductivity_outer-conductivity_inter)
%sigav is the (conductivity_inner+conductivity_outer)/2 array
%rhotoIm is the (1/conductivity_outer-1/conductivity_inter)^(-1) array
te2p=te2p-1;
nte=numel(te2p)/4;
np=numel(p)/3;np2=0;
mex_id_ = 'o int = countfaces(i int, i int[xx], i int[x], io int[x])';
[nf, np2] = BEM(mex_id_, nte, te2p, regid, np2, 4, nte, nte, 1);
t2p=zeros([3,nf]);
p2=zeros([3,np2]);
sigdiff=zeros([nf,1]);
sigav=zeros([nf,1]);
rhotoImt=sigav;
nreg=numel(unique(regid));
regid=regid-min(regid(:));
mex_id_ = 'makeBEMmesh(i int, i int[xx], i double[xx], i double[x], i int[x], io int[xx], io double[xx], io double[x], io double[x], io double[x])';
[t2p, p2, sigdiff, sigav, rhotoImt] = BEM(mex_id_, nte, te2p, p, condu, regid, t2p, p2, sigdiff, sigav, rhotoImt, 4, nte, 3, np, nreg, nte, 3, nf, 3, np2, nf, nf, nf);
t2p=t2p+1;
rhotoImt=(1./(sigav+sigdiff/2)-1./(sigav-sigdiff/2)).^-1;
rhotoIm=zeros([np2,1]);
ct=rhotoIm;
for i=1:nf
rhotoIm(t2p(:,i))=rhotoIm(t2p(:,i))+rhotoImt(i);
ct(t2p(:,i))=ct(t2p(:,i))+1;
end
rhotoIm=rhotoIm./ct;
