%final project, complex part. Alexander Wurdinger

function OLG_IDIO_AGG_RISK2

clc;
close all

global nj ny nx replrate gridx tau alpha

tic

opt_det=false;          % 1=deterministic model
opt_nosr=false;         % 1=no survival risk
opt_ny = 2;             % 1=Markov chain with number of states, ny=5,
% 2=Markov chain with ny=2 (Kr�ger-Ludwig
% calibration)
% calibration

tol = 0.1;
maxit = 20;
df = 0.3;


% guess psi 0 and 1 
%for theta=2 and rr=0 use
% psi0=[2 ;11.73];            
% psi1=[0.7; -0.4];
%for theta=2 and rr=0.6 use
% psi0=[1.99 ;1.94];            
% psi1=[0.429; 0.198];
psi0=[7.0148 ;4.6922];            %maybe kann ich das noch soft coden (initially hatte ich [1 0.9] und [alpha alpha]
psi1=[-0.1750; -0.1802];
% SOLUTION

func_calibr(opt_det,opt_nosr,opt_ny);

%VV=zeros(ny,nx,maxit);
VV=zeros(ny,nx);
for it=1:maxit,
    [psib, psir, V, gridx0] = Krusell_Smith(psi0,psi1);
    if it==1
        VV=VV+V;
     else
         %in each iteration the value functions live on different grids, so
         %we need to interpolate to be able to add them up 
         vL_intp=zeros(1,30);
         for ac=1:30
             vL_intp(ac) = interp1(gridx0(1,:),V(1,:),gridxx(1,ac),'linear');
         end
         vH_intp=zeros(1,30);
         for ac=1:30
             vH_intp(ac) = interp1(gridx0(2,:),V(2,:),gridxx(2,ac),'linear');
         end
         
         V_intp=[vL_intp; vH_intp];
         
         VV=VV+V_intp;
    end
    
    if it==1
        gridxx=gridx0;
    end
         
    
    if abs(psi0(1)-psib(1))<tol && abs(psi0(2)-psir(1))<tol && abs(psi1(1)-psib(2))<tol && abs(psi1(2)-psir(2))<tol
        display('Mambo No.5')
        VVcon=V;            %the 'real' value function if the thing actually converges
        gridxcon=gridx0;
        save vrep06t1finalcon.txt VVcon -ascii -double 
        save gridx06t1finalcon.txt gridxcon -ascii -double
        break
    end
    psi0(1)=df*psib(1)+(1-df)*psi0(1);
    psi0(2)=df*psir(1)+(1-df)*psi0(2);
    psi1(1)=df*psib(2)+(1-df)*psi1(1);
    psi1(2)=df*psir(2)+(1-df)*psi1(2);
     
end;
VV=(1/it)*VV;
psi=[psi0; psi1];
display(psi)
 save vrep06t1final.txt VV -ascii -double 
 save gridx06t1final.txt gridxx -ascii -double
end %main function 

function func_calibr(opt_det,opt_nosr,opt_ny)

global betta tetta nj jr nx ny nk pi gridy netw pens sr epsi curv pini frac pop totpop grdfac delta alpha L R replrate piz gridz gridk ETA Z nshocks 

close all

rho = 0.04;
betta = 1/(1+rho);
tetta = 1;%2;%1.1;
delta = 0.05;
alpha = 0.33;

nj=80;
jr=45;

nx=30;         % # of grid-points normal 30
curv=3.0;       % curvature of grid
grdfac=60;      % scaling factor of saving grid

% deterministic income component:
netw=1.0;
pens=0.4;%0;%
replrate=0.6;%0;%
epsi=ones(nj,1);
if (jr<nj),
    epsi(jr+1:nj)=0.0;
end;

% survival rates
if opt_nosr,
    sr = ones(nj,1);
else
    mr = readfile([],'MR.txt',3);
    sr = 1.0-mr(21:21+nj-1,1);
end;

% population and fraction living in year...
pop=zeros(nj,1);
pop(1)=100;
for jc=2:nj,
    pop(jc)=pop(jc-1)*sr(jc-1);
end;
totpop=sum(pop);

% normalize population to one:
pop=pop/totpop;
totpop=1.0;
frac=pop./totpop;

% working population
L=sum(pop(1:45));
R=1-L;

% # of income states
if (opt_det==1),
    ny = 1;
    pini = 1.0;
    gridy = 1.0;
    pi = 1.0;
else
    
    if (opt_ny==1)
        % number of income shocks
        ny = 5;
        % transition probability
        rhoeta=0.98;
        % variance of "permanent" shock
        % taken from Campbell, Viceira, Ch. 7
        vareta=0.01;
        % Markov chain:
        [pi,gridy] = markovappr(rhoeta,sqrt(vareta),2,ny);
        
        % compute invariant distribution
        pini = 1/ny*ones(ny,1);
        for tc=1:100,
            pini = pi'*pini;
        end;
        
        % take exponent and rescale income shocks such that mean is one:
        gridy=exp(gridy)';
        gridy = gridy / sum(gridy.*pini);
        
    else
        
        % Alternative -- taken from Kr�ger and Ludwig (2007) using only two
        % states
        ny = 2;
        
        % transition probability and variance
        rhoeta=0.97;
        vary=0.08;    % taken from Storesletten, Telmer, Yaron
        
        % shock
        epsil=sqrt(vary/(4.0*rhoeta*(1.0-rhoeta)));
        
        % Markov chain
        [pini,pi,gridy]=mchain(rhoeta,epsil);
    end;
    
end;

%Markov chain for aggragate risk 
piz=[0.95 0.05; 0.05 0.95];    
gridz=[1.03; 0.97];


%simulation of shocks for eta and zeta (both generated by independend
%markov chains)
nshocks=100;
zeta0=1;
% mcz = dtmc(piz);
% mce = dtmc(pi);
% Z   = simulate(mcz,nshocks,'X0',zeta0);
X=rand(nj);
eta0=zeros(1,nj);
for q=1:nj
    if X(q)>0.5
        eta0(q)=1;
    else
        eta0(q)=2;
    end
end
Z = zeros(1,nshocks);
Z(1)=zeta0;
for i=2:nshocks
    this_step_distribution = piz(Z(i-1),:);
    cumulative_distribution = cumsum(this_step_distribution);
    r = rand();
    Z(i) = find(cumulative_distribution>r,1);
end

ETA=zeros(nj,nshocks);
for q=1:nj
    e=zeros(1,nshocks);
    e(1)=eta0(q);
    
    for i=2:nshocks
    this_step_distribution = pi(e(i-1),:);
    cumulative_distribution = cumsum(this_step_distribution);
    r = rand();
    e(i) = find(cumulative_distribution>r,1);
    end
    ETA(q,:)=e;
end
    
    
    


%grid on aggreagate capital 
%for theta=2 and rr=0 use %8
%for theta=2 and rr=0.6 use %7.3108
kproj3=4;%8;%7.3108;%5.3962; 9.9680;%  %equilibrium aggregate capital taken from project 3 

kmin=0.5*kproj3;
kmax=1.5*kproj3;

nk=5;

gridk=kmin:(kmax-kmin)/(nk-1):kmax;


end     % end function func_calibr
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

function [psib, psir, V, gridx0] = Krusell_Smith(psi0,psi1)

global L R replrate delta gridx gridz tau gridk nk alpha Z

lnk1=zeros(2,nk);
for j=1:nk
    for i=1:2
        lnk1(i,j)=psi0(i)+psi1(i)*log(gridk(j));
    end
end
%display(lnk1)
k1=exp(lnk1);
%display(k1)
tau = func_pens(L,R,replrate);


% solution of household model
[gridx,gridsav,gridass,cfun,vfun] = func_hh(k1);

gass=gridx-cfun;



%simulation and regression 
[psib, psir, V, gridx0] = simulation(gridx,gridsav,gridass,cfun,vfun);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [psib, psir, V, gridx0]=simulation(gridx,gridsav,gridass,cfun,vfun)

global alpha delta ETA Z gridk gridy gridz nshocks ny nk nx

%period one
V=zeros(ny,nx);
Kstates=zeros(1,nshocks);
K=zeros(1,nshocks+1);
K(1)=gridk(3);
STK=3;
errors=zeros(1,nshocks);
for tc=1:nshocks
    wage=(1-alpha)*K(tc)^(alpha)*gridz(Z(tc));
    ret=alpha*K(tc)^(alpha-1)*gridz(Z(tc))-delta;
    gridxs=squeeze(gridx(:,:,STK,Z(tc),:));
    cfuns=squeeze(cfun(:,:,STK,Z(tc),:));
    
    % aggregation
    [Phi,PhiAss,ass,ERR] = func_aggr(gridxs,gridsav,cfuns,gridass,wage,ret);
    errors(tc)=ERR;
    %get value 
     %test=squeeze(vfun(1,:,STK,Z(tc),:));
     if tc==1
         V=V+(1/nshocks)*squeeze(vfun(1,:,STK,Z(tc),:));
     else
         %in each iteration the value functions live on different grids, so
         %we need to interpolate to be able to add them up 
         v1=squeeze(vfun(1,:,STK,Z(tc),:));
         gridx1=squeeze(gridxs(1,:,:));
         vL_intp=zeros(1,30);
         for ac=1:30
             vL_intp(ac) = interp1(gridx1(1,:),v1(1,:),gridx0(1,ac),'linear');
         end
         vH_intp=zeros(1,30);
         for ac=1:30
             vH_intp(ac) = interp1(gridx1(2,:),v1(2,:),gridx0(2,ac),'linear');
         end
         
         v_intp=[vL_intp; vH_intp];
         
         V=V+(1/nshocks)*v_intp;
         
     end  
         
     if tc==1
         gridx0=squeeze(gridxs(1,:,:));
     end
     K(tc+1)=ass;
     B = repmat(ass,[1 length(gridk)]);
     [xxx, STK]=min(abs(B-gridk));
%     K(tc+1)=gridk(STK);
      Kstates(tc)=STK;
      
      

end


sb=1;
sr=1;
KB=[];
KB1=[];
KR=[];
KR1=[];
for t=1:nshocks-1
    if Z(t)==1
        KB(sb)=K(t);
        KB1(sb)=K(t+1);
        sb=sb+1;
    else
        KR(sr)=K(t);
        KR1(sr)=K(t+1);
        sr=sr+1;
    end
end

 KB1=KB1';
 KB0=[ones(length(KB),1) KB'];
 
 KR1=KR1';
 KR0=[ones(length(KR),1) KR'];

% OLS 

psib= (KB0'*KB0)^(-1)*KB0'*KB1 ;
psir= (KR0'*KR0)^(-1)*KR0'*KR1 ;

end














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [gridx,gridsav,gridass,cfun,vfun] = func_hh(k1)

global alpha betta tetta nj nx ny nk pi gridy gridk gridz pens sr epsi curv grdfac tau replrate piz delta

disp('solution of household model');

% grids and decisions rules:
gridx = zeros(nj,ny,nk,2,nx);
gridsav = zeros(nx,1);
gridass = zeros(nj,ny,nk,2,nx);
cfun = zeros(nj,ny,nk,2,nx);
vfun = zeros(nj,ny,nk,2,nx);
vpfun = zeros(nx,1);
vptrans = zeros(nj,ny,nk,2,nx);

% savings grid: hold it constant:
maxsav=grdfac;
gridsav(2:nx)=makegrid(0.0,grdfac,nx-1,curv);
gridsav(1)=0.0;

% income states
for yc=1:ny
    for kc=1:nk 
        for zc=1:2 
            wage=(1-alpha)*gridk(kc)^(alpha)*gridz(zc);
            ret=alpha*gridk(kc)^(alpha-1)*gridz(zc)-delta;
            % cash-on-hand grid at nj:
            inc = epsi(nj)*wage*gridy(yc)*(1-tau)+(1-epsi(nj))*wage*replrate*(1-tau);
            
            % in case of no pension system, assume some minimum cash on hand:
            minx=max(inc,sqrt(eps));
            maxx=gridsav(nx)*(1.0+ret)+inc;
            gridx(nj,yc,kc,zc,:)=linspace(minx,maxx,nx);
            
            % Final period Consumption function, asset holdings, value function, including derivative
            cfun(nj,yc,kc,zc,:)=gridx(nj,yc,kc,zc,:);
            gridass(nj,yc,kc,zc,:)=(gridx(nj,yc,kc,zc,:)-inc)/(1+ret);
            vfun(nj,yc,kc,zc,:)=U(cfun(nj,yc,kc,zc,:));
            vpfun(:)=MUc(cfun(nj,yc,kc,zc,:));
            vptrans(nj,yc,kc,zc,:)=vpfun.^(-1.0/tetta);
        end
    end
end;

% Iterate Backwards
for jc=nj-1:-1:1,
    for yc=1:ny
        for kc=1:nk 
            for zc=1:2 
                 K11=k1(zc,kc);
                 
                 %today wage and ret
                 wage=(1-alpha)*gridk(kc)^(alpha)*gridz(zc);
                 ret=alpha*gridk(kc)^(alpha-1)*gridz(zc)-delta;
                 
                 %put K1 on the grid of K (bzw leave K1 aber welche stelle
                 %hat es 
                 A = repmat(K11,[1 length(gridk)]);
                 [xx, stk]=min(abs(A-gridk));
                 K1=gridk(stk);
                 
                 for xc=2:nx,
                     vp=zeros(2,2);
                     
                     for ycc=1:ny,
                         for zcc=1:2
                             % income tomorrow:
                             w1=(1-alpha)*K1^alpha*gridz(zcc);
                             ret1=alpha*K1^(alpha-1)*gridz(zcc)-delta;
                             incp1=epsi(jc+1)*w1*gridy(ycc)*(1-tau)+(1-epsi(jc+1))*w1*replrate*(1-tau);
                             
                             % Maximum cash on hand tomorrow:
                             % in case of zero savings and no pension system assume some
                             % minimum cash on hand
                             cah=max(sqrt(eps),incp1+(1.0+ret1)*gridsav(xc));
                             
                             % Interpolate derivative of value function
                             if ( cah<gridx(jc+1,ycc,stk,zcc,1)),
                                 disp('how can this be?')
                                 display(zc)
                                 display(kc)
                                 display(yc)
                                 display(ycc)
                                 display(zcc)
                                 display(ret1)
                                 display(cah)
                                 display(gridx(jc+1,ycc,stk,zcc,1))
                             end;
                             if ( cah>gridx(jc+1,ycc,stk,zcc,nx) ),
                                 % if out of bounds simply set it to decision at nx:
                                 vptr = vptrans(jc+1,ycc,stk,zcc,nx);
                             else
                                 vptr = interp1(squeeze(gridx(jc+1,ycc,stk,zcc,:)),squeeze(vptrans(jc+1,ycc,stk,zcc,:)),cah);
                             end;
                             vp(ycc,zcc)=vptr.^(-tetta);
                             %display(vp)
                         end;   %end for zcc
                     end; %end for ycc
                     
                     
                     PI=pi(yc,:)'*piz(zc,:);
                     %display(PI)
                     V=vp.*PI;
                     %display(V)
                     
                     % Euler equation: RHS
                     expvp=betta*sr(jc)*(1.0+ret)*sum(sum(V));  
                     %display(expvp)
                     % consumption
                     cfun(jc,yc,kc,zc,xc)=invut(expvp);
                     
                     % endogenous x-grid:
                     gridx(jc,yc,kc,zc,xc)=gridsav(xc)+cfun(jc,yc,kc,zc,xc);
                 end; %for xc ??
                 
                 % income (wages and pensions) in current period/age:
                 inc=epsi(jc)*wage*gridy(yc)*(1-tau)+(1-epsi(jc))*wage*replrate*(1-tau);
                 
                 % decision at minx
                 % notice: correction required for welfare calculation
                 % the above is actually slightly inefficient because xmin
                 % can be explicitly computed, then gridsav would be age and
                 % state dependent.
                 minx=max(inc,sqrt(eps));
                 if (minx<gridx(jc,yc,kc,zc,2)),
                     gridx(jc,yc,kc,zc,1)=minx;
                 else    % set it to some arbitrary fracion of x(2)
                     gridx(jc,yc,kc,zc,1)=0.9*gridx(jc,yc,kc,zc,2);
                 end;
                 
                 % Compute optimal consumption and leisure for minx
                 cfun(jc,yc,kc,zc,1)=gridx(jc,yc,kc,zc,1);
                 
                 % assets at all xc:
                 gridass(jc,yc,kc,zc,:)=(gridx(jc,yc,kc,zc,:)-inc)/(1+ret);
                 
                 % Update vfun and vpfun
                 vpfun(:)=MUc(cfun(jc,yc,kc,zc,:));
                 vptrans(jc,yc,kc,zc,:)=vpfun(:).^(-1.0/tetta);
                 
                 % Calculate value function
                 for xc=1:nx,
                     
                     v=zeros(2,2);
                     for ycc=1:ny,
                         for zcc=1:2 %fehlt
                             % income tomorrow:
                             w1=(1-alpha)*K1^alpha*gridz(zcc);
                             ret1=alpha*K1^(alpha-1)*gridz(zcc)-delta;
                             incp1=epsi(jc+1)*w1*gridy(ycc)*(1-tau)+(1-epsi(jc+1))*w1*replrate*(1-tau);
                             
                             % cah tomorrow
                             cah=max(sqrt(eps),incp1+(1.0+ret1)*gridsav(xc));
                             
                             % this should never be the case:
                             if ((cah+0.0001)<gridx(jc+1,ycc,stk,zcc,1)),
                                 warning('How can this be ?');
                             end;
                             % linear interpolation:
                             v(ycc,zcc)=func_intp(squeeze(gridx(jc+1,ycc,stk,zcc,:)),squeeze(vfun(jc+1,ycc,stk,zcc,:)),cah);
                             %display(v)
                         end %end for zcc
                     end;    % end for ycc
                     
                     PI=pi(yc,:)'*piz(zc,:);
                     %display(PI)
                     V=v.*PI;
                     %display(V)
                     % update value function
                     expv=sum(sum(V));
                     %display(expv)
                     vfun(jc,yc,kc,zc,xc)=U(cfun(jc,yc,kc,zc,xc))+betta*sr(jc)*expv;
                 end;    % end for xc
                 
            end %end for zc
            
        end %end for kc
        
    end;    % end for yc
    
end;    % end for jc


% ---------------------------------------------------------------------
    function fv = func_intp(x,func,xp)
        
        
        n = length(x);
        if ( xp>x(n) ),
            % fv = func(n);
            fv=func_extrapol(x(n-1),x(n),func(n-1),func(n),xp);
        elseif (xp<x(1)),
            % fv = func(1);
            fv=func_extrapol(x(1),x(2),func(1),func(2),xp);
        else
            fv = interp1(x,func,xp);
        end;
        
    end
% ---------------------------------------------------------------------


% ---------------------------------------------------------------------
    function y=func_extrapol(x1,x2,y1,y2,x)
        
        % simple linear extrapolation
        
        m = (y2-y1)/(x2-x1);
        y = y1 + m*(x-x1);
        
    end
% ---------------------------------------------------------------------

end     % end function func_hh
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [Phi,PhiAss,ass,ERR]=func_aggr(gridxs,gridsav,cfuns,gridass,wage,ret)

global  nj nx ny pi gridy pens sr epsi pini frac totpop replrate tau

%disp('aggregation and cross-sectional measure');

% Compute Cross sectional distributions and aggregate variables
Phi = zeros(nj,ny,nx);          % distribution of assets conditional by age and shock
PhiAss = zeros(nx,1);           % distribution of assets

% Distribution of newborns over cash at hand
for yc=1:ny
    
    % income (wages and pensions) in current period/age:
    inc=epsi(1)*wage*gridy(yc)*(1-tau)+(1-epsi(1))*wage*replrate*(1-tau);
    
    % initial cash-on-hand:
    cahini=inc;
    
    [vals,inds]=basefun(gridxs(1,yc,:),cahini,nx);
    Phi(1,yc,inds(1))=vals(1)*pini(yc)*frac(1);
    Phi(1,yc,inds(2))=vals(2)*pini(yc)*frac(1);
end;

for jc=2:nj
    TT = zeros(ny,nx,ny,nx);    % transfer function
    
    for xc=1:nx
        for yc=1:ny
            for ycc=1:ny
                
                % income (wages and pensions) in current period/age:
                inc=epsi(jc)*wage*gridy(ycc)*(1-tau)+(1-epsi(jc))*wage*replrate*(1-tau);
                
                % cash on hand: x=a*(1+r)+y = s(-1)*(1+r)+y;
                cah=inc+(1.0+ret)*gridsav(xc);
                
                [vals,inds]=basefun(gridxs(jc,ycc,:),cah,nx);
                
                TT(ycc,inds(1),yc,xc)=vals(1)*pi(yc,ycc);
                TT(ycc,inds(2),yc,xc)=vals(2)*pi(yc,ycc);
            end;    
        end;    
    end;    
    
    for xc=1:nx
        for yc=1:ny
            for xcc=1:nx
                for ycc=1:ny
                    % transfer distribution:
                    Phi(jc,ycc,xcc)=Phi(jc,ycc,xcc)+Phi(jc-1,yc,xc)*TT(ycc,xcc,yc,xc)*sr(jc-1);
                end;
            end;
        end;
    end;
    
end;    % end for jc


% Check that for each country distribution sums to 1
sumprob=sum(sum(sum(Phi(:,:,:))));
ERR=0;
if ( ( sumprob < 0.999 ) || ( sumprob > 1.001) ),
    ERR=1-sumprob;
    %beep; beep; beep;
    %warning('distribution fucked up');
end;

% Check if Grid is Big enough
sumprob=sum(sum(Phi(:,:,nx)));
if (sumprob > 0.001 ),
    beep; beep; beep;
    warning('grid too small -- increase your grid');

end;

ass=0.0;
cons=0.0;


% aggregation
for jc=1:nj
    for yc=1:ny
        for xc=1:nx,
            PhiAss(xc)=PhiAss(xc)+Phi(jc,yc,xc);
            
            % asset holdings = capital stock in general equilibrium
            ass=ass+totpop*Phi(jc,yc,xc)*gridsav(xc);
            
            cons=cons+totpop*Phi(jc,yc,xc)*cfuns(jc,yc,xc);
            
            %whatever does do, if i need them. we will see
%             lab=lab+totpop*Phi(jc,yc,xc)*gridy(yc)*epsi(jc);
%             ret=ret+totpop*Phi(jc,yc,xc)*gridy(yc)*(1.0-epsi(jc));
%             
        end;
    end;
end;


% ---------------------------------------------------------------------
    function [vals,inds]=basefun(grid_x,x,nx)
        % this subroutine returns the values and the indices of the two basis
        % functions that are positive on a given x in the grid_x
        
        % MF function to lookup the current position
        i=lookuppp(grid_x,x,0);
        
        if ( (i+1)>nx),
            inds(1)=nx;
            inds(2)=nx;
            vals(2)=0.0;
            vals(1)=1.0;
        elseif (i==0),
            inds(1)=1;
            inds(2)=1;
            vals(1)=1.0;
            vals(2)=0.0;
        else
            inds(1)=i;
            inds(2)=i+1;
            dist = grid_x(i+1)-grid_x(i);
            vals(2)=( x-grid_x(i) )/dist;
            vals(1)=( grid_x(i+1)-x )/dist;
        end;
        
    end 	% end function basefun
% ---------------------------------------------------------------------




end     % end function func_aggr
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

function [mpk,Y] = func_mpk(ass, L)

global alpha

Y = ass.^alpha * L.^(1-alpha);
ky = ass./Y;
mpk = alpha * ky.^(-1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pini,pi,gridy]=mchain(rhoeta,epsil)

% Transition Probabilities
pi=rhoeta*ones(2,2);
pi(1,2)=1.0-rhoeta;
pi(2,1)=1.0-rhoeta;

% Initial Distribution
pini=0.5*ones(2,1);

gridy=zeros(2,1);
gridy(1)=exp(1.0-epsil);
gridy(2)=exp(1.0+epsil);
gridy=2.0*gridy/(sum(gridy));

end  % end function mchain
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function u = U(c)
global tetta

% utility
if (abs(tetta-1-0)<sqrt(eps)),
    u = log(c);
else
    u = c.^(1.0-tetta)/(1.0-tetta);
end;
end     % end function U
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function muc=MUc(c)
global tetta

% maringal utility
muc = c.^(-tetta);
end     % end function MUc
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function invut=invut(marg)
global tetta

% invert utility for c
invut=marg.^(-1.0/tetta);
end     % end function invut
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function grd = makegrid(x1,x2,n,c)
% makes curved grid according to curvature parameter c
scale=x2-x1;
grd(1)=x1;
grd(n)=x2;
for i=2:n-1
    grd(i)=x1+scale*((i-1.0)/(n-1.0))^c;
end;
end     % end function makegrid
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

function tau = func_pens(L,R,replrate)

tau = replrate*R ./ (L + replrate * R);

end
% ----------------------------------------------

function [v0] = value0(Phi,gridass,cfun,vfun)

global  nj nx ny  frac betta

disp('life-cycle profiles')

asslife=zeros(nj,1);
conslife=zeros(nj,1);


for jc=1:nj,
    for yc=1:ny
        for xc=1:nx,
            %asslife(jc)=asslife(jc)+Phi(jc,yc,xc)*gridass(jc,yc,xc)/frac(jc);
            conslife(jc)=conslife(jc)+Phi(jc,yc,xc)*cfun(jc,yc,xc)/frac(jc);
            
        end;
    end;
end;

v0=0;
for jc=1:nj
   v0=v0+betta^(jc-1)*conslife(jc);
end
display(v0)
end     % end function lcprofile
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++