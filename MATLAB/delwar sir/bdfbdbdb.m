clear;clc;
x=0:.5:4;
h=x(2)-x(1);
y=[0,13,33,39.5,40,40,36,15,0]./60;
oddsum=0;evensum=0;
for n=1:1:length(x)
    if mod(n,2)==0
        evensum=evensum+y(n)*2;
    else
        oddsum=oddsum+y(n)*4;
    end
end
sum=(h/3)*(oddsum+evensum);
fprintf('Distance between two stations %.3f',sum);