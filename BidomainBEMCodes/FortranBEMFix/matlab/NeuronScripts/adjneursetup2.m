function BEM=adjneursetup2(tri,p,sigdiff,sigav,nlen,fmmquad,c,conv,ind_p_mem)

nt=numel(tri)/3;
np=numel(p)/3;
%%%Define integral equations
%Calculate the n x grad N for the Hypersingular int operator
BEM.condudiff=cat(2,sigdiff(:),sigdiff(:),sigdiff(:))';
BEM.condudiff=sparse(1:3*nt,1:3*nt,BEM.condudiff(:));
BEM.conduav=cat(2,sigav(:),sigav(:),sigav(:))';
nv=1:nt;
whattriag=cat(2,nv(:),nv(:),nv(:))';
BEM.conduav=sparse(1:3*nt,1:3*nt,BEM.conduav(:));
BEM.NttoN=sparse(1:3*nt,tri(:)',ones([1 3*nt]));
ncgradn=getncrossgradN(p',tri');
BEM.ncrossgrad{1}=sparse(whattriag(:)',1:3*nt,reshape(ncgradn(:,1,:),[1 3*nt]));
BEM.ncrossgrad{2}=sparse(whattriag(:)',1:3*nt,reshape(ncgradn(:,2,:),[1 3*nt]));
BEM.ncrossgrad{3}=sparse(whattriag(:)',1:3*nt,reshape(ncgradn(:,3,:),[1 3*nt]));
%Get the near-field for the Green's function kernal
[BEM.AnearS0,BEM.AnearD0,BEM.AnearS,BEM.AnearD,col,row,BEM.area]=nearmat2(tri,p,nlen);

ID=grammsurf(tri',p',sigav); % sigav int Ni Nj dS
BEM.AnearD=-ID+BEM.NttoN'*BEM.condudiff*BEM.AnearD'*BEM.NttoN; %Complete A matrix near field
%Finish the calculation of the hyper singular operator
BEM.AnearS0=...
BEM.NttoN'*BEM.condudiff*(...
                  BEM.ncrossgrad{1}'*BEM.AnearS0*BEM.ncrossgrad{1}...
                 +BEM.ncrossgrad{2}'*BEM.AnearS0*BEM.ncrossgrad{2}...
                 +BEM.ncrossgrad{3}'*BEM.AnearS0*BEM.ncrossgrad{3})*BEM.NttoN;


BEM.Amat=BEM.AnearD;
%Add A matrix and Nmem to get the total left-hand side
BEM.Amat(:,ind_p_mem)=BEM.AnearD(:,ind_p_mem)+ c*BEM.AnearS0(:,ind_p_mem)*conv(ind_p_mem,ind_p_mem);

%Prepare the preconditioning for our matrix inversion later
[~,~,BEM.Ach,BEM.Alin,BEM.rfmm,BEM.nhatfmm]=precor(tri,p,1,1,fmmquad);

nnz(BEM.AnearD)
clear D0p S0p;
BEM.nfmm=numel(BEM.rfmm)/3;

BEM.preconeul=sparse(1:np,1:np,1./diag(BEM.AnearD));
BEM.preconAmat=sparse(1:np,1:np,1./diag(BEM.Amat));

BEM.Amat=BEM.preconAmat*BEM.Amat;
BEM.AnearD=BEM.preconeul*BEM.AnearD;
end

function lap=grammsurf(t2p,p,sigmaav)
np=numel(p)/3;
nt=numel(t2p)/3;
col=ones([nt*9 1]);
row=ones([nt*9 1]);
val=zeros([nt*9 1]);
for i=1:nt
col(9*(i-1)+1:9*i)=t2p(i,[1 1 1  2 2 2  3 3 3 ]);
row(9*(i-1)+1:9*i)=t2p(i,[1 2 3  1 2 3  1 2 3 ]);
BEM.area=0.04166666666666666*norm(cross(...
p(t2p(i,2),:)-p(t2p(i,1),:),p(t2p(i,3),:)-p(t2p(i,1),:)));

val(9*(i-1)+1:9*i)=sigmaav(i)*[2,1,1,1,2,1,1,1,2]*BEM.area;
end
lap=sparse(col,row,val,np,np);
end

