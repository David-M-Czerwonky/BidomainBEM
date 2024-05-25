function Eprimary=Eprimfmm(errfmm,rs,js,robs)
ns=numel(rs)/3;
nobs=numel(robs)/3;
errfmm=10^-12;
mex_id_ = 'computeEprimary(i double[x], i double[xx], i double[xx], i int[x], i double[xx], i int[x], o double[xx])';
[Eprimary] = BEM(mex_id_, errfmm, rs, js, ns, robs, nobs, 1, 3, ns, 3, ns, 1, 3, nobs, 1, 3, nobs);
end






