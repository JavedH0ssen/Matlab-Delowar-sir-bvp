clc;clear;
f=@(x) 1/(4*x+5);
a=0;b=5;n=20;
h=(b-a)/n;oddsum=0;evensum=0;sum3=0;sum=0;sum6=0;
for i=1:n-1
    s=f(a+i*h);
    sum=sum+s;
    if mod(i,2)==0
        evensum=evensum+s;
    else
        oddsum=oddsum+s;
    end
    if mod(i,3)==0
        sum3=sum3+s;
    end
    if mod(i,6)==0
        sum6=sum6+s;
    end
end
nonthird=sum-sum3;non6_3=sum3-sum6;non3odd=oddsum-non6_3;
simp1_3=h/3*(f(a)+f(b)+4*oddsum+2*evensum);
weddle=(3*h/10)*(f(a)+f(b)+evensum+sum6+5*non3odd+6*non6_3);
fprintf("result of integration by\nsimpsons 1/3=%5f \n" + ...
    "weddle =%.5f ",simp1_3,weddle);