%% Q3 
% Work in progress. There is a problem with my cheb approximation. Using a
% plonomial interpolation (and cheb nodes), the value function converges to
% the same function as before ! 
clear;clc;
%% Step 1, discretizes the state space
tic
kss=50;
kmin=10; 
kmax=kss;
n=199;
%K=kmin:(kmax-kmin)/n:kmax; 
%p=size(K,2);
theta=0.679;
h=1;
delta=0.013;
beta=0.988;

%%
chevn1=zeros(1,n);

for j=1:n
    b=(2*j-1)/(2*n);
   chevn1(j)=cos(b*pi);
end

chevn=flip(chevn1);

cb=(chevn+1)*((kmax-kmin)/2)+kmin;

%%
x=1:11;
K=cb;
p=size(K,2);
V=zeros(p);
V1=ones(1,p);

%initial guess V=0;

m=zeros(p);
M=zeros(p);

for i=1:p
    for j=1:p
        m(i,j)=K(i)^(1-theta)*h^theta+(1-delta)*K(i)-K(j);
    end
end

% Step 4, unfeasible k'
for i=1:p
    for j=1:p
        if m(i,j)>0
            M(i,j)=log(m(i,j));
        else 
            M(i,j)=-1000;
        end
    end
end
        
% Step 5 Matrix chi and updated value function 

X=M+beta*V;


for i=1:p
    V1(i)=max(X(i,:));
end

%%
d=5;
pp5=polyfit(K,V1,5);

f=@(x) pp5(6)+pp5(5)*x+pp5(4)*x.^2+pp5(3)*x.^3+pp5(2)*x.^4+pp5(1)*x.^5;

V=zeros(1,p);
%%
eps=0.1;
s=1;
while double(any(abs(V1-V)>eps)) > 0

    V=V1;
    
m=zeros(p);
M=zeros(p);

for i=1:p
    for j=1:p
        m(i,j)=K(i)^(1-theta)*h^theta+(1-delta)*K(i)-K(j);
    end
end

%
for i=1:p
    for j=1:p
        if m(i,j)>0
            X(i,j)=log(m(i,j))+beta*(pp5(6)+pp5(5)*K(j)+pp5(4)*K(j)^2+pp5(3)*K(j)^3+pp5(2)*K(j)^4+pp5(1)*K(j)^5);
        else 
            X(i,j)=-1000;
        end
    end
end

for i=1:p
    V1(i)=max(X(i,:));
end


d=5;
pp5=polyfit(K,V1,5);



s=s+1;

    if s==1000
        break
    end

end
toc

plot(K,V1)
