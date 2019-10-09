%% Quetion 2, HW3
clear;clc;
%parameters
kapA=[5 5];
vA=[1 1];
sigA=[0.8 0.8];
nlA=[5.5 3.5];
nhA=[0.5 2.5];
ZA=[1 1];
thetaA=[0.6 0.6];
KA=[2 2];
klA=[1 1];
khA=[1 1];
lambA=[0.95 0.84];
phiA=[0.2 0.2]; 

%%
Cl=zeros(1,2);
Ch=zeros(1,2);
Hl=zeros(1,2);
Hh=zeros(1,2);
Wage=zeros(1,2);
Rate=zeros(1,2);
C=zeros(1,2);
y=zeros(1,2);
for i=1:2
kap=kapA(i);     
v=vA(i);         
sig=sigA(i);     
nl=nhA(i);       %swap high and low bc of missspecification in the model
nh=nlA(i);       
Z=ZA(i);         
theta=thetaA(i); 
K=KA(i);         
lamb=lambA(i);   
phi=phiA(i);     
kl=klA(i);
kh=khA(i);



%exogenous variables, cl,ch,hl,hh


syms cl ch hl hh w r

%Euler equations of high and low typ
EE1 = (1-phi)*cl^(-sig)*lamb*(w*hl*nl)^(-phi)*w*nl-kap*hl^(1/v)==0;
EE2 = (1-phi)*ch^(-sig)*lamb*(w*hh*nh)^(-phi)*w*nh-kap*hh^(1/v)==0;

%prices, already use market clearing. Goods marke clears by walras law
R = (1-theta)*Z*K^(-theta)*(nl*hl+nh*hh)^(theta)==r; 
W = theta*Z*K^(1-theta)*(nl*hl+nh*hh)^(theta-1)==w;

%BCs
BCL = lamb*(w*hl*nl)^(1-phi)+r*kl==cl;
BCH = lamb*(w*hh*nh)^(1-phi)+r*kh==ch;

%numerical solve system of equations
[s1,s2, s3, s4, s5, s6] = vpasolve([EE1,EE2,R,W,BCL,BCH],[ cl, ch, hl, hh, w, r]);


Cl(i)=double(s1);
Ch(i)=double(s2);
Hl(i)=double(s3);
Hh(i)=double(s4);
Wage(i)=double(s5);
Rate(i)=double(s6);

y(i)=Z*K^(1-theta)*(nl*Hl(i)+nh*Hh(i))^theta;   % check if goods market clears
C(i)=Cl(i)+Ch(i);
end



