%HW5 
clear;clc;
%% Q1,2
gam=0.6; %0.8
n=10000000;
mu=[1 1];
sig=[1 0;0 1]; %[1 0;0 1]; %[1 0.5; 0.5 1]; %[1 -0.5; -0.5 1];

zk = mvnrnd(mu,sig,n);

lnk=zk(:,1)';
lnz=zk(:,2)';

%%
k=exp(lnk);
z=exp(lnz);

%% 2

y=z.*k.^gam;

%% 3

K=sum(k);

zz=zeros(1,n);
ke=zeros(1,n);

for i=1:n
    zz(i)=(z(1)/z(i))^(1/(gam-1));
end

ke(1)=K/sum(zz);

for i=2:n
    ke(i)=ke(1)*zz(i);
end


%% 4 

gap=k-ke;

%% 5

Ya=sum(y);

ye=z.*ke.^gam;

Ye=sum(ye);


OGR=(Ye/Ya-1)*100;


%% Q3
% 3 do it 1000 times
p=1000;
sOGR=zeros(1,p);
for j=1:p
% 1 
ss=10000; %100; %1000; %10000 %100000              %samplesize
szk=datasample(zk,ss,1);
slnk=szk(:,1)';
slnz=szk(:,2)';


cvm=cov(slnk,slnz);     %covariance matrix 


%% 2 

sk=exp(slnk);
sz=exp(slnz);

%% 2.2

sy=sz.*sk.^gam;

%% 2.3

sK=sum(sk);

szz=zeros(1,ss);
ske=zeros(1,ss);

for i=1:ss
    szz(i)=(sz(1)/sz(i))^(1/(gam-1));
end

ske(1)=sK/sum(szz);

for i=2:ss
    ske(i)=ske(1)*szz(i);
end


%% 2.4 

sgap=sk-ske;

%% 2.5

sYa=sum(sy);

sye=sz.*ske.^gam;

sYe=sum(sye);


sOGR(j)=(sYe/sYa-1)*100;

end


%% descriptive statistics

histogram(sOGR)

M=median(sOGR);

%% 4

lb=0.9*OGR;
ub=1.1*OGR;


in= sOGR(sOGR>lb & sOGR<ub);

prob=size(in,2)/p;


