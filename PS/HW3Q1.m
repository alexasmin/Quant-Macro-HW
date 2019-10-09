%% HW3 quant macro 

clear;clc;      
h=0.31;
theta=0.67;
delta=0.0625;      %define variables with the analytically obtained 
beta=0.98;         %steady state values 
z0=1.6298;
z1=z0*2;
y0=1;
k0=4;
kss=8;
css=1.5;
c0=0.75; 
%%
%first period
k=k0;

%guess c1
count2=1;
c=1;
while css-c>0.02  %first while loop updates the guess for ct to get closer to the 
k=k0;             %correct solution. As the second loop stops when k is close
seqk=[];          %enough to the new steady state k, the outer loop of updated
seqc=[];          %guess will stop when c is close enough at css at the same time
seqi=[];
seqy=[];
    
    
c=1+0.0005*count2;  %increase guess in each iteration 
    

count=1;

while kss-k>0.1   %inner loop, calculates the transition until k is close enough
                  %to the new steady state value
y=k^(1-theta)*((z1*h)^theta);
k1=y-c+(1-delta)*k;

%Euler Equation
c1=beta*c*((1-delta)+(1-theta)*(k1^(-theta))*(z1*h)^theta);
seqk(count)=k;
seqc(count)=c;
seqi(count)=y-c;   %saving the sequence of all important variables
seqy(count)=y;
k=k1;
c=c1;
count=count+1;
end
count2=count2+1;

end
%%
seqk(count)=k;
seqc(count)=c;
y=k^(1-theta)*((z1*h)^theta);
seqi(count)=y-c;
seqy(count)=y;
count2=count2+1;
size(seqk)
%%
plot(seqk)
title('Transition Path of Capital')
%%
plot(seqc)
hold on 
plot(seqi)
hold on 
plot(seqy)
legend('Consumption',...
       'Investment',...
       'Output',  'location','best')              
title('Transitoin Paths')

hold off


%% d %%%%%%%%%%

clear; 
h=0.31;
theta=0.67;
delta=0.0625;
beta=0.98;
z0=1.6298;
z1=z0*2;
y0=1;
k0=4;
kss=4;      %old staedy state will be reached again 
css=0.75;
c0=0.75; 
%%
%first period
k=k0;

%guess c1 the same as in c 

c=1.0275; 

seqk=[];
seqc=[];
seqi=[];
seqy=[];
count=1;

%%
for i=1:10;                          %calculate the transition for the first
    y=k^(1-theta)*((z1*h)^theta);    %ten periods, same transition as in c
    k1=y-c+(1-delta)*k;
    %Euler Equation

    c1=beta*c*((1-delta)+(1-theta)*(k1^(-theta))*(z1*h)^theta);
    seqk(i)=k;
    seqc(i)=c;
    seqi(i)=y-c;
    seqy(i)=y;
    k=k1;
    c=c1;
    count=count+1;

end

%%
count2=1;
c=1;
                     
%%                
                  %now the same as in c only for the fall after the shocj
                  %in 10. IMPORTANT: Stopping criterions reversed
while c-css>0.02  %first while loop updates the guess for ct to get closer to the 
k=seqk(10);             %correct solution. As the second loop stops when k is close
seqk(11:end)=[];          %enough to the new steady state k, the outer loop of updated
seqc(11:end)=[];          %guess will stop when c is close enough at css at the same time
seqi(11:end)=[];
seqy(11:end)=[];
    
    
c=1-0.0005*count2;  %increase guess in each iteration 
    

count=1;
while k-kss>0.1
    
y=k^(1-theta)*((z0*h)^theta);
k1=y-c+(1-delta)*k;
%Euler Equation

c1=beta*c*((1-delta)+(1-theta)*(k1^(-theta))*(z0*h)^theta);
seqk(10+count)=k;
seqc(10+count)=c; 
seqi(10+count)=y-c;
seqy(10+count)=y;
k=k1;
c=c1;
count=count+1;
end
count2=count2+1;
end

%%
seqk(10+count)=k;
seqc(10+count)=c;
y=k^(1-theta)*((z0*h)^theta);
seqi(10+count)=y-c;
seqy(10+count)=y;

%%
plot(seqk)
title('Transition Path of Capital')
%%
plot(seqc)
hold on 
plot(seqi)
hold on 
plot(seqy)
legend('Consumption',...
       'Investment',...
       'Output',  'location','best')              
title('Transition Paths')
hold off





