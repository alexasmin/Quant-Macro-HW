%% Q3 
% Work in progress. Each iteration takes a very long time and the value
% function does not seem to converge
clear;clc;
%% Step 1, discretizes the state space
kss=50;
kmin=10; 
kmax=kss;
n=200;
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
d=3;
thetas=ones(d+1,1);

thetas(1)=sum(V1)/p;

%%
for j=1:d;

thetas(j+1)=cheb(p, j, K, V1);

end

V=zeros(1,p);
%%
s=1;
while double(any(abs(V1-V(1,:))>eps)) > 0

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
            X(i,j)=log(m(i,j))+beta*(thetas(1)+thetas(2)*K(j)+thetas(3)*(2*K(j)^2-1)+thetas(4)*(4*K(j)^3-3*K(j)));
        else 
            X(i,j)=-1000;
        end
    end
end

for i=1:p
    V1(i)=max(X(i,:));
end


d=3;
thetas=ones(d+1,1);

thetas(1)=sum(V1)/p;


for j=1:d;

thetas(j+1)=cheb(p, j, K, V1);

end
s=s+1;
    if s==10
        break
    end
end