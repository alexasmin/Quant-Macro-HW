%% Q2 b 
clear;clc;
%parameter values
kap=5;
v=1;
sig=0.8;
nl=5.5;
nh=0.5;
Z=1;
theta=0.6;
KA=2;
KLA=1;
KHA=1;
lamb=0.95;
phi=0.2; 
kapB=5;
vB=1;
sigB=0.8;
nlb=3.5;
nhb=2.5;
thetaB=0.6;
KB=2;
KLB=1;
KHB=1;
lambB=0.84;
phiB=0.2; 

%exogenous variables
syms cla cha hla hha wa ra kla kha clb chb hlb hhb wb rb klb khb

%country A
%Euler equations 
EE1A = (1-phi)*cla^(-sig)*lamb*(wa*hla*nl)^(-phi)*wa*nl-kap*hla^(1/v)==0;
EE2A = (1-phi)*cha^(-sig)*lamb*(wa*hha*nh)^(-phi)*wa*nh-kap*hha^(1/v)==0;
EK1A = nl*ra*kla^(nl-1)-rb==0;
EK2A = nh*ra*kha^(nh-1)-rb==0;

%prices, already use market claearing. Goods marke claears by walras law
RA = (1-theta)*Z*(kla+kha+KB-klb-khb)^(-theta)*(nl*hla+nh*hha)^(theta)==ra; 
WA = theta*Z*(kla+kha+KB-klb-khb)^(1-theta)*(nl*hla+nh*hha)^(theta-1)==wa;

%BCs
BCLA = lamb*(wa*hla*nl)^(1-phi)+ra*kla^nl+rb*(KLA-kla)==cla;
BCHA = lamb*(wa*hha*nh)^(1-phi)+ra*kha^nh+rb*(KHA-kha)==cha;

%Country B
%Euler equations
EE1B = (1-phiB)*clb^(-sigB)*lambB*(wb*hlb*nlb)^(-phiB)*wb*nlb-kapB*hlb^(1/vB)==0;
EE2B = (1-phiB)*chb^(-sigB)*lambB*(wb*hhb*nhb)^(-phiB)*wb*nhb-kapB*hhb^(1/vB)==0;
EK1B = nlb*ra*klb^(nlb-1)-rb==0;
EK2B = nhb*ra*khb^(nhb-1)-rb==0;

%prices, already use market clbearing. Goods marke clbears by wbalras lawb
RB = (1-thetaB)*Z*(klb+khb+KB-kla-kha)^(-thetaB)*(nlb*hlb+nhb*hhb)^(thetaB)==rb; 
WB = thetaB*Z*(klb+khb+KB-kla-kha)^(1-thetaB)*(nlb*hlb+nhb*hhb)^(thetaB-1)==wb;

%BCs

BCLB = lambB*(wb*hlb*nlb)^(1-phiB)+rb*klb^nlb+ra*(KLB-klb)==clb;
BCHB = lambB*(wb*hhb*nhb)^(1-phiB)+rb*khb^nlb+ra*(KHB-khb)==chb;


%% numerical solve the system
[s1,s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16] = vpasolve([EE1A,EE2A,EK1A,EK2A,EE1B,EE2B,EK1B,EK2B, RA,WA,RB,WB,BCLA,BCHA BCLB, BCHB],[ cla, cha, hla, hha, wa, ra, kla, kha, clb, chb, hlb, hhb, wb, rb, klb, khb]);

%%
CLA=double(s1);
CHA=double(s2);
HLA=double(s3);
HHA=double(s4);
WA=double(s5);
RA=double(s6);
KLA=double(s7);
KHA=double(s8);
CLB=double(s9);
CHB=double(s10);
HLB=double(s11);
HHB=double(s12);
WB=double(s13);
RB=double(s14);
KLB=double(s15);
KHB=double(s16);


%% check market clearing
y=Z*(KLA+KHA+2-KLB-KHB)^(1-theta)*(nl*HLA+nh*HHA)^theta;

C=CHA+CLA;


