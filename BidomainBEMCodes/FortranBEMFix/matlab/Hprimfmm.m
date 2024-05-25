function Hout=Hprimfmm(iprecEp,rs,js,robs)
nsource=numel(rs)/3;
ntarget=numel(robs)/3;
Hout=zeros([3 ntarget]);
mex_id_ = 'Hprim(i int, i int, i double[xx], i double[xx], i double[xx], io double[xx], i double)';
[Hout] = BEM(mex_id_, nsource, ntarget, rs, js, robs, Hout, iprecEp, 3, nsource, 3, nsource, 3, ntarget, 3, ntarget);

end





