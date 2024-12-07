f=@(x,y) (x+y+x*y);
x0=0;
y0=1;
h=0.025;
Xn=0.2;
y=y0;
for x=x0:h:Xn
    if x==.1 ||x==.2
        fprintf('y(%0.1f)=%0.4f\n',x,y);
    end
    y=y+h*f(x,y);
 
end