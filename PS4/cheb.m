function [ theta ] = cheb( nodes, degree, x, y)
%chebishev theta given number nodes, the degree and the x and y values

j=nodes;
d=degree;

zahlerele=ones(j,1);

for i=1:j;

 zahlerele(i)= y(i)*chebyshevT(d, x(i));
end
zahler=sum(zahlerele); 

nennerele=ones(j,1);

for i=1:j;
    
    nennerele(i)= chebyshevT(d, x(i))*chebyshevT(d, x(i));
    
end

nenner=sum(nennerele);

theta=zahler/nenner;
end


