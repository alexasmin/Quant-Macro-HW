function vpp1 = evalvp(cons,x)

global betta tetta r g gridx vpfun epsi probepsi ne nx 

vpp1 = zeros(ne,1);
lam1=zeros(ne,1);
for ec=1:ne,
    xp1 = (x-cons)*(1+r)/(1+g)+epsi(ec);
    lam=vpfun.^(-1/tetta);
    lam1(ec) = (func_intp(gridx,lam,xp1));
    vpp1=lam1.^(-tetta);
end;
vpp1 = sum(vpp1.*probepsi);

end
