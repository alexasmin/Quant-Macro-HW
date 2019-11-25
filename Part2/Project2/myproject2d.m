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
ct = zeros(nt,1);
xt = zeros(nt,1);
at = zeros(nt,1);
yt = zeros(nt,1);
et = zeros(nt,1);
eulert = zeros(nt,1);

% random numbers:
rand('seed',0);
probet=rand(nt,1);
indet=ceil(probet*ne);

outx = 0;
for tc=1:nt,
    yt(tc) = epsi(indet(tc));
    xt(tc) = at(tc)+yt(tc);
    ct(tc) = func_intp(gridx,cons,xt(tc));
    chkoutx=false;
    if (ct(tc)==max(cons)),
        outx=outx+1;
        chkoutx=true;
    end;
    if (tc<nt),
        at(tc+1)=(xt(tc)-ct(tc))*(1+r)/(1+g);
    end;
    
    % error evaluation
    if (ct(tc)<xt(tc) && chkoutx==false),
        margu=margutil(ct(tc));
        et(tc)=abs(foc(ct(tc),xt(tc))/margu);
    end;
end;
% -------------------------------------------------------------------------

avget=mean(et);

fracoutx=outx/nt;
%if ( routx>0.01 ),
%    beep; beep; beep;
%    warning('grid too small, enlarge your grid');
    disp(['fraction of points outside grid is ', num2str(fracoutx)]);
%end;

nnt=nt-dt+1;
disp(['maximum Euler equation error is : ', num2str(max(et(dt+1:nt)))]);
disp(['mean Euler equation error is : ', num2str(mean(et(dt+1:nt)))]);
I = (et==0);
rbc = sum(I(dt+1:nt))/nnt;
disp(['borrowing constraint binds in ', num2str(rbc), ' cases']);

figure;
plot([dt+1:nt],ct(dt+1:nt),'LineWidth',2);
xlabel('time');
ylabel('c_t');
title('consumption over time');

figure;
plot([dt+1:nt],xt(dt+1:nt),'LineWidth',2);
xlabel('time');
ylabel('x_t');
title('cash-on-hand over time');

figure;
plot([dt+1:nt],et(dt+1:nt),'LineWidth',2);
xlabel('time');
ylabel('e_t');
title('euler error over time');



