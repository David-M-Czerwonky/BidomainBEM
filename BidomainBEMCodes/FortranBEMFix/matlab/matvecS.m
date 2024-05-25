function y=matvecS(Achlin,rfmm,nfmm,Anear,x)
ytemp=Achlin*(x/(4*pi));
errfmm=10^-12;
mex_id_ = 'fmmS(i double[x], i double[x], i double[xx], i int[x], o double[x])';
[ytemp] = BEM(mex_id_, errfmm, ytemp, rfmm, nfmm, 1, nfmm, 3, nfmm, 1, nfmm);
y=Achlin'*ytemp+Anear*x;
end


