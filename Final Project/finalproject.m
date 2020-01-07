% Final Project, ALexander Wurdinger
clc;clear;

%% 1.2 

%set up
alpha=0.3;
beta=0.99^40;
tau=0;
lam=0.5;
T=50000;
b=0;
maxiter=100;
tol=0.0001;
ww=0.2;

%shocks 

lnzeta=normrnd(0,0.13,1,T);
%lnrohh=normrnd(0,0.5,1,T);

lnz=zeros(1,T);
%lnroh=zeros(1,T);

for t=1:T
    if lnzeta(t)<=0
        lnz(t)=-0.13;
    else
        lnz(t)=0.13;
    end
    
%     if lnrohh(t)<=1
%         lnroh(t)=-0.5;
%     else
%         lnroh(t)=0.5;
%     end
end

z=exp(lnz);
% income states zeta 
lnzz=[-0.13; 0.13];
zz=exp(lnzz);


%% income shocks and roh

mat=load('etafinal.txt');
lneta=mat(:,2);
probeta=mat(:,3);
eta=exp(lneta);

lnroh=[-0.5; 0.5];  % both states assumed to have equal probability
roh=exp(lnroh);
%% savings rate (and Phi) 

phis=zeros(2,11);


for i=1:2
   for j=1:11
    phis(i,j)=1/(1+((1-alpha)/(alpha*(1+lam)*roh(i)))*(lam*eta(j)+tau*(1+lam*(1-eta(j)))));
   end
   
end
 

phi=sum(phis(1,:)*probeta*0.5)+sum(phis(2,:)*probeta*0.5);   


s=beta*phi/(1+beta*phi);

%% simulation of k 

%steady state of k in logs 
lnkss=(log(s)+log(1-tau)+log(1-alpha)-log(1+lam))/(1-alpha);
kss=exp(lnkss);
%simulation with productivity shocks 
lnk=zeros(1,T);
lnk(1)=lnkss;

%%
for t=2:T
    lnk(t)=log(s)+log(1-tau)+log(1-alpha)-log(1+lam)+log(z(t))+alpha*lnk(t-1);
end


%% Krusell Smith 

kmin=log(0.5)+lnkss;
kmax=log(1.5)+lnkss;

n=5;

gridk=kmin:(kmax-kmin)/(n-1):kmax;

% guess psi 0 and 1 
psi0=[-2.205 ;-1.1945];  %maybe kann ich das noch soft coden 
psi1=[alpha; alpha];

%% loop 
for q=1:maxiter
% final period 

lnk1=zeros(2,n);
for j=1:n
    for i=1:2
        lnk1(i,j)=psi0(i)+psi1(i)*gridk(j);
    end
end

k1=exp(lnk1);

a=zeros(2,n);
c1=zeros(2,n);
c2=zeros(2,n);
s=zeros(2,n);
W=0;
%%
for j=1:n
    for i=1:2 %today prod shock
        %R=alpha*k1(i,j)^(alpha-1);
        w=(1-alpha)*exp(gridk(j))^(alpha)*zz(i);
        
        WW=zeros(1,44);
        BB=zeros(1,4);
        d=1;
        dd=1;
        for m=1:2            %tomorrow prod shock
            for p=1:2        %tomorrow ret shock
                for e=1:11   %tomorrow income shock
                    %wt1=(1-alpha)*k1(i,j)
                    %W=W+(0.5*0.5*probeta(e))*(eta(e)*(1-alpha)*k1(i,j)*zz(m))/(alpha*k1(i,j)^(alpha-1)*zz(m)*roh(p));
                    wt1=(1-alpha)*k1(i,j)^(alpha)*zz(m);
                    etaR=eta(e)/(alpha*k1(i,j)^(alpha-1)*zz(m)*roh(p));
                    WW(d)=wt1*etaR*0.5*0.5*probeta(e);
                    W=sum(WW);
                    d=d+1;
                end
                RR=alpha*k1(i,j)^(alpha-1)*zz(m)*roh(p);
                BB(dd)=(wt1/RR)*0.5*0.5;
                B=sum(BB);
            end
        end
        
        a(i,j)=((1-tau)*w-(1+lam)*(tau/beta)*B-(1/beta)*lam*(1-tau)*W)/(1+(1/beta));
        if a(i,j)<0
            a(i,j)=0.00001;    %non negativity. Important. Otherwise log(s) can be complex !!
        end
        s(i,j)=a(i,j)/((1-tau)*w);
        c1(i,j)=(1-tau)*w-a(i,j);
        %c2(i,j)=a(i,j)*R+lam*b; %fehlt natürlich noch zeug
    end
end

%% simulation 
lnksim=zeros(1,T);

lnksim(1)=gridk(3);   %start in steady state
stk=3;
for t=2:T
    if lnz(t)==0.13
        st=2;
    else
        st=1;
    end
%     stk=find(gridk == lnksim(t-1));
%     lnksim1=log(s(st,stk))+log(1-tau)+log(1-alpha)-log(1+lam)+log(z(t))+alpha*lnksim(t-1);
%     A = repmat(lnksim1,[1 length(gridk)]);
%    [xx, ind]=min(abs(A-gridk));
%    lnksim(t)=gridk(ind);

%or 

    
    lnksim1=log(s(st,stk))+log(1-tau)+log(1-alpha)-log(1+lam)+log(z(t))+alpha*lnksim(t-1);
    %interpolating k to get the right state for next period s
    A = repmat(lnksim1,[1 length(gridk)]);
   [xx, stk]=min(abs(A-gridk));
   lnksim(t)=lnksim1;

end


%% regression

sb=1;
sr=1;
KB=[];
KB1=[];
KR=[];
KR1=[];
for t=501:T-1
    if lnz(t)==0.13
        KB(sb)=lnksim(t);
        KB1(sb)=lnksim(t+1);
        sb=sb+1;
    else
        KR(sr)=lnksim(t);
        KR1(sr)=lnksim(t+1);
        sr=sr+1;
    end
end
%%
% KB1=KB(2:end)';
 KB1=KB1';
 KB0=[ones(length(KB),1) KB'];
% 
% KR1=KR(2:end)';
 KR1=KR1';
 KR0=[ones(length(KR),1) KR'];

% OLS 

psib= (KB0'*KB0)^(-1)*KB0'*KB1 ;
psir= (KR0'*KR0)^(-1)*KR0'*KR1 ;

if abs(psi0(1)-psir(1))<tol && abs(psi0(2)-psib(1))<tol && abs(psi1(1)-psir(2))<tol && abs(psi1(2)-psib(2))<tol
    display('Mambo No.5')
    break
end


psi0(1)=ww*psir(1)+(1-ww)*psi0(1);
psi0(2)=ww*psib(1)+(1-ww)*psi0(2);
psi1(1)=ww*psir(2)+(1-ww)*psi1(1);
psi1(2)=ww*psib(2)+(1-ww)*psi1(2);


end

%% EU 

%simulate also a sequence of shocks for roh and for eta 

zsim=z;

lnrohhsim=normrnd(0,0.5,1,T);

lnrohsim=zeros(1,T);

for t=1:T

    
    if lnrohhsim(t)<=1
        lnrohsim(t)=-0.5;
    else
        lnrohsim(t)=0.5;
    end
end

rohsim=exp(lnrohsim);


et=rand(1,T);
etasim=zeros(1,T);

bounds=zeros(1,12);

% sequence of etas

for i=2:12
    bounds(i)=bounds(i-1)+probeta(i-1);
end

for t=1:T
    for j=1:11
        if et(t)>=bounds(j) 
            etasim(t)=eta(j);
        end
    end
end



%% ex post
% ksim=exp(lnksim);
% c1sim=zeros(1,T-1);
% c2sim=zeros(1,T-1);
% EU=zeros(1,T-1);
% for t=1:T-1
%     R=alpha*ksim(t)^(alpha-1)*rohsim(t)*zsim(t);
%     w=(1-alpha)*ksim(t)^alpha*rohsim(t);
%     b=tau*w*((1-lam)/(1+lam));
%     c1sim(t)=(1-tau)*w-ksim(t+1);
%     c2sim(t)=ksim(t)*R+lam*etasim(t)*w*(1-tau)+(1-lam)*b;
% end
% 
% %%
% for t=1:T-1
%     EU(t)=log(c1sim(t))+log(c2sim(t));
% end
% 
% 
% avgEU=1/(T-1-500)*sum(EU(501:end));


%% ex ante second try
ksim=exp(lnksim);
c1sim=zeros(1,T-1);
Ec2sim=zeros(44,1);
V=zeros(1,T-1);
EU=zeros(1,T-1);
Beta=beta/(1+beta);
prob=repmat(probeta,[4,1])*0.5*0.5;
for t=1:T-1
    w=(1-alpha)*ksim(t)^alpha*rohsim(t);
    c1sim(t)=(1-tau)*w-ksim(t+1);
    
    d=1;
    for m=1:2            %tomorrow prod shock
        for p=1:2        %tomorrow ret shock
            for e=1:11   %tomorrow income shock
                
                R=alpha*ksim(t+1)^(alpha-1)*roh(p)*zz(m);
                w1=(1-alpha)*ksim(t+1)^alpha*zz(m);
                b=tau*w1*((1-lam)/(1+lam));
                Ec2sim(d)=ksim(t)*R+lam*eta(e)*w1*(1-tau)+(1-lam)*b;
                d=d+1;
            end
        end
    end
    
    util=log(Ec2sim);
    V(t)=sum(util.*prob);
    
                    
end

%%
for t=1:T-1
    EU(t)=(1-Beta)*log(c1sim(t))+Beta*V(t);
end


avgEU=1/(T-1-500)*sum(EU(501:end));

%% consumption variation
V1=-1.7658;
V2=-1.88;

g=exp((V2-V1)/beta)-1;










