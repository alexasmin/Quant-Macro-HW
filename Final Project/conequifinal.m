%Plot average utilities and calculate consumption equivalent

tetta=2;
rho = 0.04;
betta = 1/(1+rho);

 V1=load('vrep06t1finalcon.txt','-ascii');
 V0=load('vrep0t1finalcon.txt','-ascii');
% 
 X1=load('gridx06t1finalcon.txt','-ascii');
 X0=load('gridx0t1finalcon.txt','-ascii');

plot(X1(2,:),V1(2,:))
plot(X0(2,:),V0(2,:)) 
 
%%
v0L_intp=zeros(1,30);
    for ac=1:30
        v0L_intp(ac) = interp1(X0(1,:),V0(1,:),X1(1,ac),'linear');
    end
v0H_intp=zeros(1,30);
    for ac=1:30
        v0H_intp(ac) = interp1(X0(2,:),V0(2,:),X1(2,ac),'linear');
    end
    
v0_intp=[v0L_intp; v0H_intp];

%%
%g=(V1./v0_intp).^(1/(1-tetta))-1;

g=exp((V1-v0_intp)./betta)-1;
%%
plot(X1(2,1:20),g(2,1:20))
