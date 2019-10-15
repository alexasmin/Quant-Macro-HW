%% HW4 
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
s=0;               %counting the number of iterations
%% Step 6 checking if fixed point is reached 

while double(any(abs(V1-V(1,:))>eps)) > 0                           %any(abs(V1-V(1,:))>eps)

    if s>0;
        V=repmat(V1,p,1);
    end
    
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
    V1(i)=max(X(i,:));
end

s=s+1;
%if double(any(abs(V1-V(1,:))>eps)) == 0
    %b=1;
    %break
end
toc   
%end
%%

plot(K,V1)

