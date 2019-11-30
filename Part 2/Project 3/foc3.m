function fval = foc3(cons,x,t,yc)

global betta tetta r g gridx vpfun epsi probepsi ne nx sr

vpp1 = evalvp3(cons,x,t,yc);
margu = MUc(cons);
fval = margu - betta*sr(t)*(1.0+r) * vpp1;

end