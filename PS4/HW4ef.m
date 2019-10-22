%policy function iterations
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
s=0;          %counting the number of iterations
ss=0;
done=0;
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

sp=[50 100 150]; %starting points for the policy function iteration
st=[5 10 20 50]; %number of policy iterations w/o changing the policy function


for b=1:sp(1) %run the value function iterations normally in the beginning
     if s>0;
        V=repmat(V1,p,1);
    end        
% Step 5 Matrix chi and updated value function 

X=M+beta*V;


for i=1:p
    [V1(i), g(i)]=max(X(i,:));
end

s=s+1;

end




%% 
%no iterate using the policy function. Use different amounts of iterations
%before updating the policy function 
while double(any(abs(V1-V(1,:))>eps)) > 0
    V=repmat(V1,p,1);
    
for b=1:st(4)
    
    for i=1:p
        V1(i)=log(K(i)^(1-theta)*h^theta+(1-delta)*K(i)-K(g(i)))+beta*V(1,g(i));
    end
    ss=ss+1;
    if double(any(abs(V1-V(1,:))>eps)) == 0
        done=1;
        break
    end
    V=repmat(V1,p,1);
    
end

if done==1;
    break
end
s=s+1;

X=M+beta*V;


for i=1:p
    [V1(i), g(i)]=max(X(i,:));
end

s=s+1;

end
toc
%%
plot(K,V1)