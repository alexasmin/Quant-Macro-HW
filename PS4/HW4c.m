%% c concavity

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
done=1;   %1 as long as one value is not close enough 
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
        
% Step 5.1 Matrix chi and updated value function 
%%
Xj=0;
Xj1=0;
for i=1:p                           %instead of calculating the whole X matrix
    for j=1:p-1                     %and then max each row, calculate each 
        Xj=M(i,j)+beta*V(i,j);      %entry and stop if Xj>Xj+1 and take Xj
        Xj1=M(i,j+1)+beta*V(i,j+1); %as the max. 
        if Xj>Xj1                   %Note: for this problem and in matlab it seems
            V1(i)=Xj;               %to actually slow down the problem, as matrix
            break                   %addition is not costly. Maybe for more 
        end                         %grid entries the loop safes time over the 
                                    %max function 
    end
end
        






s=s+1;
%if double(any(abs(V1-V(1,:))>eps)) == 0
    %b=1;
    %break
   
end
toc
%%
plot(V1)