% Consumption equivalent 

tetta=2;

 v0=load('v0_2.txt','-ascii');
% vT=load('vT_2.txt','-ascii');
% 
 x0=load('x0_2.txt','-ascii');
% xT=load('xT_2.txt','-ascii');

% v0=load('v0_1.txt','-ascii');
% vT=load('vT_1.txt','-ascii');
% 
% x0=load('x0_1.txt','-ascii');
% xT=load('xT_1.txt','-ascii');

vTdec=load('vT_dec.txt','-ascii');
xTdec=load('xT_dec.txt','-ascii');

v0L_intp=zeros(1,30);
    for ac=1:30
        v0L_intp(ac) = interp1(x0(1,:),v0(1,:),xTdec(1,ac),'linear');
    end
v0H_intp=zeros(1,30);
    for ac=1:30
        v0H_intp(ac) = interp1(x0(2,:),v0(2,:),xTdec(2,ac),'linear');
    end
    
v0_intp=[v0L_intp; v0H_intp];

g=(vTdec./v0_intp).^(1/(1-tetta))-1; %tetta=/0

%g=exp(vT./v0_intp)-1; %tetta=1, log case noooo

plot(xTdec(2,:),g(2,:))


