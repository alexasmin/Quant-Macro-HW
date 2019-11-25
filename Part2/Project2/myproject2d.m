% Problem d 
%% Fintite time 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% project 2
% solution of household problem for T = \infty
%function project2

clear all
close all

global betta tetta r g gridx vpfun epsi probepsi ne nx min_cons

% -------------------------------------------------------------------------
% SETTINGS
maxit = 100; 
tol=1e-4;
nt=1100;    % periods for simulation
dt=100;     % periods discarded
min_cons=1.0e-08;

% parameters
r = 0.02;
rho = 0.03;
g = 0.01;
tetta = 1;
betta = 1/(1+rho);
T=60;

% grid
nx=50;              % # of grid-points
curv=3.0;           % curvature of grid
xmax = 30;          % scaling factor of saving grid
xmin = sqrt(eps);
gridx=makegrid(xmin,xmax,nx,curv);
gridx=gridx';

% income shocks
ne = 7;
varepsi = 0.01;
muepsi = -varepsi/2;
% [epsi,probepsi] = qnwnorm(ne,muepsi,varepsi);
% mat=[[1:ne]',epsi,probepsi];
% save rn.txt mat -ascii -double -tabs
mat=load('rn.txt');
epsi=mat(:,2);
probepsi=mat(:,3);
epsi=exp(epsi);
if (abs(sum(epsi.*probepsi)-1.0)>sqrt(eps)),
    error('random numbers fucked up');
end;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%% SOLUTION 
    cons=zeros(nx,T);
    cons(:,T)=gridx;
   
        %%
     for t=T-1:-1:1
       vpfun = margutil(cons(:,t+1));
       for xc=1:nx,
           % check binding constraint:
           mincons=gridx(xc);
           mu = foc(mincons,gridx(xc));
           if (mu>=0.0),
               cons(xc,t)=mincons;
           else
               [cons(xc,t),fval] = fzero('foc',cons(xc,t),[],gridx(xc));
               cons(xc,t)=max(cons(xc,t),min_cons);
           end;
       end;    
     end
    % update vpfun
   
% -------------------------------------------------------------------------
%%
figure;
plot(gridx,cons(:,1),'LineWidth',2);
xlabel('x');
ylabel('c');
title('consumption policy');
%%
%[gridx,cons]



% -------------------------------------------------------------------------
% SIMULATION
ct = zeros(T,1);
xt = zeros(T,1);
at = zeros(T,1);
yt = zeros(T,1);
et = zeros(T,1);
eulert = zeros(T,1);

% random numbers:
rand('seed',0);
probet=rand(T,1);
indet=ceil(probet*ne);

outx = 0;
for tc=1:T,
    yt(tc) = epsi(indet(tc));
    xt(tc) = at(tc)+yt(tc);
    ct(tc) = func_intp(gridx,cons(:,tc)',xt(tc));
    chkoutx=false;
    if (ct(tc)==max(cons(:,tc))),
        outx=outx+1;
        chkoutx=true;
    end;
    if (tc<T),
        at(tc+1)=(xt(tc)-ct(tc))*(1+r)/(1+g);
    end;
    
    % error evaluation
    if (ct(tc)<xt(tc) && chkoutx==false),
        margu=margutil(ct(tc));
        et(tc)=abs(foc(ct(tc),xt(tc))/margu);
    end;
end;
%% -------------------------------------------------------------------------
figure;
plot(ct);
xlabel('time');
ylabel('c_t');
title('consumption over time');



