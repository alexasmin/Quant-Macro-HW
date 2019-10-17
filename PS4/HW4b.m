%% b taking into account monotonicity

clear;clc;

%% Step 1, discretizes the state space
tic
kss=50;
kmin=10; 
kmax=kss;
n=199;
K=kmin:(kmax-kmin)/n:kmax; 
p=size(K,2);
theta=0.679;
h=1;
delta=0.013;
beta=0.988;
%% Step 2, Initial guess for V
eps=ones(1,p)*0.1;  %toleranz value
V=zeros(p);
V1=ones(1,p);
g=ones(1,p);
s=0;
%% 
    
% Step 3, Return matrix

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
    [V1(i), g(i)]=max(X(i,:));
end

s=s+1;
%if double(any(abs(V1-V(1,:))>eps)) == 0
    %b=1;
    %break
    
%% Step 6     
    

    
% Step 3, Return matrix

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
        
% Step 5.1 Matrix chi and updated value function %%%%%%%%%%%%%%%
while double(any(abs(V1-V(1,:))>eps)) > 0                           
    
    if s>0;
        V=repmat(V1,p,1);
    end

%%
gg=zeros(1,p);       %all the g,gg,ggg are to save the actual position of the max
ggg=zeros(1,p);
for i=1:p
    Xr=zeros(1,p-g(i)+1);     %to avoid calculating the whole X matrix, only
    t=1;
    for j=g(i):p              %calculate for position gresater or equal then 
    Xr(t)=M(i,j)+beta*V(i,j); %the policy function of last period tells you
    t=t+1;
    end
    [V1(i), gg(i)]=max(Xr);    %then maximize each row and safe
    ggg(i)=g(i)-1+gg(i);; %calculate the true position of the arg max
end

g=ggg;


s=s+1;
%if double(any(abs(V1-V(1,:))>eps)) == 0
    %b=1;
    %break
    
end
toc
%%   
plot(V1)


