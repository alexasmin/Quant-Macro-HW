%Project 2 

%% Problem c 
clear;clc;

maxit = 100; 
tol=1e-4;
nt=1100;    % periods for simulation
dt=100;     % periods discarded
min_cons=1.0e-08;

% parameters
r = 0.02;
rho = 0.03;
gr = 0.01;
tetta = 1;
betta = 1/(1+rho);

% grid
nx=50;              % # of grid-points
curv=3.0;           % curvature of grid
xmax = 30;          % scaling factor of saving grid
xmin = sqrt(eps);
gridx=makegrid(xmin,xmax,nx,curv);
x=gridx';
pp=size(x,1);

% income shocks
ne = 7;
varepsi = 0.01;
muepsi = -varepsi/2;
% [epsi,probepsi] = qnwnorm(ne,muepsi,varepsi);
% mat=[[1:ne]',epsi,probepsi];
% save rn.txt mat -ascii -double -tabs
mat=load('rn.txt');
epsi=mat(:,2);
probepsi=mat(:,3);
epsi=exp(epsi);
if (abs(sum(epsi.*probepsi)-1.0)>sqrt(eps)),
    error('random numbers fucked up');
end;

meaneps=mean(epsi.*probepsi);

%% no borrowing in the state space
a=zeros(42,1);
aaa=x-meaneps;
aa=aaa(10:end);    %das hier noch hard coden bitten 
a(2:end)=aa;
%% with borrowing n the state space
% a=x-meaneps;

%%

p=size(a,1);



%% fit epsilons with the grid 
% epsifit=zeros(1,ne);
% epspo=zeros(1,ne);
% mid=zeros(1,nx-1);
% for i=1:nx-1
%     mid(i)=(x(i)+x(i+1))/2;
% end
% 
% 
% for i=1:nx-1
%     for j=1:ne
%         if (epsi(j)>=x(i)) && (epsi(j)<=x(i+1))
%             if epsi(j)>=mid(i)
%                 epsifit(j)=x(i+1);
%                 epspo(j)=i+1;
%             else 
%                 epsifit(j)=x(i);
%                 epspo(j)=i;
%             end
%         end
%     end
% end


% Solution 

%% Step 2, Initial guess for V
V=zeros(p,p,ne);
w=zeros(1,p);
V1=ones(1,p,ne);
VV=zeros(1,ne);
g=ones(1,p,ne);
s=0;               %counting the number of iterations
%%     
% Step 3, Return matrix

e=mean(epsi.*probepsi);   % expected shock 

m=zeros(p,p,ne);
M=zeros(p,p,ne);

    for i=1:p
        for j=1:p
            for t=1:ne
%                 if x(j)-epsi(t)>=0
                   m(i,j,t)=a(i)+epsi(t)-a(j)/(1+r);
%                 else
%                     m(i,j,t)=-1;
%                 end
            end
        end
    end

if tetta==1
        
    % Step 4, unfeasible x' that is if the borrowing constraint is binding
    % in this set up this includs cases of negative consumption whcich have
    % to be exluded anyway
    st=0;
    for i=1:p
        for j=1:p
            for t=1:ne
                if m(i,j,t)>0
                    M(i,j,t)=log(m(i,j,t));
                    st=st+1;
                else
                    M(i,j,t)=-1000;
                end
            end
        end
    end
        
else
    
    % Step 4, unfeasible x'
    for i=1:p
        for j=1:p
            for t=1:ne
                if m(i,j,t)>0
                    M(i,j,t)=(1/(1-tetta))*m(i,j,t).^(1-tetta);
                else
                    M(i,j,t)=-1000;
                end
            end
        end
    end

end

MM1=M(:,:,1);

% Loop
%% Step 5
for it=1:maxit
    disp(['iteration # ', num2str(it)]);
   if it>1;
       V=V1; 
   end
       for i=1:p
           for j=1:ne
               VV(j)=V1(1,i,j);
           end
           w(i)=mean(VV.*probepsi');
       end
      
      W=repmat(w,[p,1,ne]);
  


X=zeros(p,p,ne);
for t=1:ne
    X(:,:,t)=M(:,:,t)+(betta*(1+gr)^(1-tetta))*W(:,:,t);
end

%first find the max of each row in the first layer and safe them in the
%first ro of V1K, then the same for q layers and safe them in q rows. That
%is find the best combination of x and x' for each epsilon

for t=1:ne
    for i=1:p
       [V1(1,i,t), g(1,i,t)]=max(X(i,:,t));
    end
end

%%
   
if double(any(abs(V1(1,:,:)-V(1,:,:))>tol)) == 0
    break;
end
  

end

if (it==maxit),
    warning('increase # of iters for specified tolerance');
end;


% %% consumption policy function
% c=zeros(p*ne,2);
% 
% l=1;
% for j=1:ne
%    for i=1:p
%         c(l,1)=a(i)+epsi(j);
%         c(l,2)=a(i)+epsi(j)-a(j)/(1+r);
%         l=l+1;
%     end
% end
% 
% c=c';
% 
% [temp, order] = sort(c(1,:));
% cc = c(:,order);

%%
plot(a,V1(:,:,1))
plot(a,g(:,:,1))



%% Howard improvment

V=zeros(p,p,ne);
w=zeros(1,p);
V1=ones(1,p,ne);
VV=zeros(1,ne);
g=ones(1,p,ne);

sp=[1 100 150]; %starting points for the policy function iteration
st=[5 10 20 50]; %number of policy iterations w/o changing the policy function

%%

for b=1:sp(1) %run the value function iterations normally in the beginning

   
       for i=1:p
           for j=1:ne
               VV(j)=V1(1,i,j);
           end
           w(i)=mean(VV.*probepsi');
       end
      
      W=repmat(w,[p,1,ne]);
  


X=zeros(p,p,ne);
for t=1:ne
    X(:,:,t)=M(:,:,t)+(betta*(1+gr)^(1-tetta))*W(:,:,t);
end

%first find the max of each row in the first layer and safe them in the
%first ro of V1K, then the same for q layers and safe them in q rows. That
%is find the best combination of x and x' for each epsilon

for t=1:ne
    for i=1:p
       [V1(1,i,t), g(1,i,t)]=max(X(i,:,t));
    end
end

end

%%
for it=1:maxit
    
    V=V1;
    
     for i=1:p
           for j=1:ne
               VV(j)=V1(1,i,j);
           end
           w(i)=mean(VV.*probepsi');
      end
      
      W=repmat(w,[p,1,ne]);
      
      %policy function iteration
      for b=1:st(4)
          
         
          for i=1:p
              for j=1:ne
                  V1(1,i,j)=M(i,g(1,i,j),j)+(betta*(1+gr)^(1-tetta))*W(1,g(1,i,j),t);
              end
          end
          
          if double(any(abs(V1(1,:,:)-V(1,:,:))>tol)) == 0  %break process if convergence is reached
             done=1;                                         %while iterating ploicy function
             break
          end  
          V=V1;
      end
      %end of policy function iteration 
    
     V=V1;
    
     %re calculate new policy function 
     
      for i=1:p
           for j=1:ne
               VV(j)=V1(1,i,j);
           end
           w(i)=mean(VV.*probepsi');
       end
      
      W=repmat(w,[p,1,ne]);
  


X=zeros(p,p,ne);
for t=1:ne
    X(:,:,t)=M(:,:,t)+(betta*(1+gr)^(1-tetta))*W(:,:,t);
end

%first find the max of each row in the first layer and safe them in the
%first ro of V1K, then the same for q layers and safe them in q rows. That
%is find the best combination of x and x' for each epsilon

for t=1:ne
    for i=1:p
       [V1(1,i,t), g(1,i,t)]=max(X(i,:,t));
    end
end
    
           if double(any(abs(V1(1,:,:)-V(1,:,:))>tol)) == 0  %break process if convergence is reached
             done=1;                                         %while iterating ploicy function
             break
           end  
    
    
end
    
    
    
    
    
    
