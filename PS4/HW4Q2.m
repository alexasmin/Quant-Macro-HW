%% Labour supply 
clear;clc;

%% Step 1, discretizes the state space
kss=50;
kmin=10; 
kmax=kss;
n=199;
K=kmin:(kmax-kmin)/n:kmax; 
p=size(K,2);
hmin=0.01;
hmax=1;
nn=99;
H=hmin:(hmax-hmin)/nn:hmax;
q=size(H,2);
theta=0.679;
delta=0.013;
beta=0.988;
kap=5.24;
v=2;

%%
eps=ones(1,p)*0.1;  %toleranz value
V=zeros(p,p,q);
V1=ones(1,p);
% Step 3, Return matrix

m=zeros(p,p,q);
M=zeros(p,p,q);



% creating 3D "Matrix" that maps all possible combinations of k,k',h


for i=1:p
    for j=1:p
        for t=1:q
          m(i,j,t)=K(i)^(1-theta)*H(t)^theta+(1-delta)*K(i)-K(j); 
        end
    end
end

%%  
% again getting rid of values for k' and h that lead to negative
% consumption then adding disutility of work to the return matrix
%stays constant, so no need to put it in the loop
for i=1:p
    for j=1:p
        for t=1:q
            if m(i,j,t)>0
                M(i,j,t)=log(m(i,j,t))-kap*(H(t)^(1+(1/v))/(1+(1/v)));
            else 
                M(i,j,t)=-1000;
            end
        end
    end
end

% no loop over all iterations of the chi matrix
s=0; 
while double(any(abs(V1-V(1,:,1))>eps)) > 0  %stopping createrion

    if s>0;
        V=repmat(V1,[p,1,q]);
    end
%%
%q layers of chi matrixes
X=zeros(p,p,q);
for t=1:q
    X(:,:,t)=M(:,:,t)+beta*V(:,:,q);
end
    
%%
%first find the max of each row in the first layer and safe them in the
%first ro of V1K, then the same for q layers and safe them in q rows. That
%is find the best combination of k and k' for each h
V1k=zeros(q,p);
for t=1:q
    for i=1:p
       V1k(t,i)=max(X(i,:,t));
    end
end

%%
%Now find the max of each coloumn of V1k. That is find the best value of h.
%Yields the new Values function 
for j=1:p
    V1(j)=max(V1k(:,j));
end
    
s=s+1;    
end    

%%
plot(K,V1)