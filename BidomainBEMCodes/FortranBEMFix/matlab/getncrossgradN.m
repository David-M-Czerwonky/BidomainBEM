function ncgradn=getncrossgradN(p,tri);
nt=numel(tri)/3;
np=numel(p)/3;
nhat=zeros(nt,3);
ncgradn=zeros(3,3,nt);
iv=eye(3);
iv(end+1,:)=0;
coeff=ones([4 4]);
coeff(end,end)=0;
for i=1:nt
    q=p(tri(i,:),:);
dir(1,:)=q(2,:)-q(1,:);
dir(2,:)=q(3,:)-q(1,:);
nhat(i,:)=reshape([dir(1,2)*dir(2,3)-dir(1,3)*dir(2,2) ...
                   dir(1,3)*dir(2,1)-dir(1,1)*dir(2,3) ...
                   dir(1,1)*dir(2,2)-dir(1,2)*dir(2,1)],[1 3]);
nhat(i,:)=nhat(i,:)/norm(nhat(i,:));
coeff(1:3,1:end-1)=q;
coeff(4,1:end-1)=nhat(i,:);
abcd=coeff\iv;
ncgradn(:,:,i)=abcd(1:3,:)';

for j=1:3
ncgradn(j,:,i)=reshape([nhat(i,2)*ncgradn(j,3,i)-nhat(i,3)*ncgradn(j,2,i) ...
                        nhat(i,3)*ncgradn(j,1,i)-nhat(i,1)*ncgradn(j,3,i) ...
                        nhat(i,1)*ncgradn(j,2,i)-nhat(i,2)*ncgradn(j,1,i)],[1 3]);
end
end


