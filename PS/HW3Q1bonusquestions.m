%% BONUS questions on question 1

clear;clc;  
tauc=0.5;
h=0.31;
theta=0.67;
delta=0.0625;      %define variables with the analytically obtained 
beta=0.98;         %steady state values 
z0=1.6298;
z1=z0*2;
y0=1;
k0=4;
kss=8;
css=1.5/(1+tauc);
c0=0.75/(1+tauc); 



%%
%first period
k=k0;

%guess c1

c=0.685; %0.934

seqk=[];
seqc=[];
seqi=[];
seqy=[];
count=1;
%%
while kss-k>0.1
    
y=k^(1-theta)*((z1*h)^theta);
k1=y-(1+tauc)*c+(1-delta)*k;
%Euler Equation

c1=beta*c*((1-delta)+(1-theta)*(k1^(-theta))*(z1*h)^theta);
seqk(count)=k;
seqc(count)=c;
seqi(count)=y-c;
seqy(count)=y;
k=k1;
c=c1;
count=count+1;
end
%%
seqk(count)=k;
seqc(count)=c;
y=k^(1-theta)*((z1*h)^theta);
seqi(count)=y-c;
seqy(count)=y;
size(seqk)
%%
plot(seqk)


