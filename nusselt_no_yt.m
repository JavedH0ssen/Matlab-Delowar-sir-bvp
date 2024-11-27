function   nusselt_no_yt
clear all
clc
close all
c       =  {'-b',  '--b', '.:b'};
c1       =  {'-r',  '--r', '.:r'};
Pr       =  [6.8   6.8   6.8];
M        =  [0  0.3  0.5];
Kp       =  [0.   0.1  0.2   0.3 0.4  0.5];
A        =  [0.2    0.2    0.2];
epsi     =  [0.5  0.5    0.5 ];
Nb       =  [0.1    0.1    0.1 ];
Nt       =  [0.1    0.1     0.1];
Rd       =  [0.2     0.2   0.2 ];
Le       =  [0.5    0.5    0.5];
Lb       =  [0.5    0.5    0.5];
s        =  [0.1  0.1  0.1 ];
Ec       =  [0    0    0];
Pe       = [0.5    0.5    0.5];
Kr       = [0.2   0.2   0.2];
sigma1   =[0.1    0.1   0.1];
    function ysol = bvpex1(x,y)        %bvpex1(ind variable, dependent variable)

        yy1    = -y(1)*y(3) +  y(2)^2 - 1 - A(i) + A(i)*((x/2)*y(3)+y(2))....
            +(M(i)+Kp(j))*(y(2)-1);

        yy2     = -(1/(1+(4/3)*Rd(i)))*(Pr(i)*y(1)*y(5) + Nb(i)*y(5)*y(7)+....
            Nt(i)*y(5)*y(5) +Pr(i)*(Ec(i)*y(3)*y(3) + s(i)*y(4)-(x/2)*A(i)*y(5)));

        yy3    =  -(Nt(i)/Nb(i))*yy2 - Le(i)*Pr(i)*y(1)*y(7) + (x/2)*Le(i)*Pr(i)*A(i)*y(7)...
            +Le(i)*Pr(i)*Kr(i)*y(6);

        yy4    = -Lb(i)*Pr(i)*y(1)*y(9) + Pe(i)*((y(8)+sigma1(i))*yy3 + y(7)*y(9))+....
            (x/2)*Lb(i)*Pr(i)*A(i)*y(9);

        ysol   = [y(2); y(3); yy1; y(5); yy2 ; y(7) ; yy3 ; y(9) ; yy4];
    end

% Here I define residual of boundary conditions
    function res = bcex1(y0, yinf)

        res  =  [y0(1)  ; y0(2) - epsi(i) ; yinf(2)-1 ;  y0(4)-1  ; yinf(4)  ;  y0(6)-1 ; yinf(6) ; y0(8)-1 ; yinf(8)];
        %y0(1)= f(1),y0(2)=f'(1),yinf(2)=f'(infinity)

    end
fdp=zeros(6,1)
for i=1:3
    for j=1:6

        % Initial values
        sol1 = bvpinit(linspace(0,10,25), [1 1 0 0 0 0 0   0   0]);
        % sol1 = bvpinit(linspace(xa,xb,n), [initiall guesses for f f' f'' etc]);

        % bvp4c takes two functions defined below and give solution in structure
        % form
        sol = bvp5c(@bvpex1, @bcex1, sol1);

        % x values
        x = sol.x;%creates mesh

        % solution , (y, y', y'', theta, theta', )
        y = sol.y;

        % find all y(0), y'(0), y''(0) so on
        value = deval(sol, 0);  %evaluate(f,f',f'', at t=1)


        % plot any solution, I plot y''
        fdp(j)=y(3,1);
    end
    hold on
    figure (1)
    plot(Kp, fdp,c{i}, 'linewidth', 2.5);

    xlabel('Kp')
    ylabel('f''"(0)')
    title('Pr = 6.8, Rd=Kr=0.2, Nb=Nt=s=0.1, Le=Lb=Pe=0.5, Ec=0','fontsize',10,'fontweight','bold','color','black')
    hold on


end




A        =  [0.5    0.5    0.5];


    function ysol1 = bvpex2(x,y)        %bvpex1(ind variable, dependent variable)

        yy1    = -y(1)*y(3) +  y(2)^2 - 1 - A(i) + A(i)*((x/2)*y(3)+y(2))....
            +(M(i)+Kp(j))*(y(2)-1);

        yy2     = -(1/(1+(4/3)*Rd(i)))*(Pr(i)*y(1)*y(5) + Nb(i)*y(5)*y(7)+....
            Nt(i)*y(5)*y(5) +Pr(i)*(Ec(i)*y(3)*y(3) + s(i)*y(4)-(1/2)*A(i)*x*y(5)));

        yy3    =  -(Nt(i)/Nb(i))*yy2 - Le(i)*Pr(i)*y(1)*y(7) + (1/2)*Le(i)*Pr(i)*A(i)*x*y(7)...
            +Le(i)*Pr(i)*Kr(i)*y(6);

        yy4    = -Lb(i)*Pr(i)*y(1)*y(9) + Pe(i)*((y(8)+sigma1(i))*yy3 + y(7)*y(9))+....
            (1/2)*Lb(i)*Pr(i)*x*A(i)*y(9);

        ysol1   = [y(2); y(3); yy1; y(5); yy2 ; y(7) ; yy3 ; y(9) ; yy4];
    end

% Here I define residual of boundary conditions
    function res1 = bcex2(y0, yinf)

        res1  =  [y0(1)  ; y0(2) - epsi(i) ; yinf(2)-1 ;  y0(4)-1  ; yinf(4)  ;  y0(6)-1 ; yinf(6) ; y0(8)-1 ; yinf(8)];
        %y0(1)= f(1),y0(2)=f'(1),yinf(2)=f'(infinity)

    end

for i=1:3
    for j=1:6
        % Initial values
        sol1 = bvpinit(linspace(0,3,25), [1 1 0 0 0 0 0  0   0]);
        % sol1 = bvpinit(linspace(xa,xb,n), [initiall guesses for f f' f'' etc]);

        % bvp4c takes two functions defined below and give solution in structure
        % form
        sol = bvp5c(@bvpex2, @bcex2, sol1);

        % x values
        x = sol.x;%creates mesh

        % solution , (y, y', y'', theta, theta', )
        y = sol.y;

        % find all y(0), y'(0), y''(0) so on
        value = deval(sol, 0);  %evaluate(f,f',f'', at t=1)


        fdp(j)=y(3,1);
    end
    hold on
    figure (1)
    plot(Kp, fdp,c1{i}, 'linewidth', 2.5);

    xlabel('Kp')
    ylabel('f(0)')
    title('Pr = 6.8, Rd=Kr=0.2, Nb=Nt=s=0.1, Le=Lb=Pe=0.5, Ec=0','fontsize',10,'fontweight','bold','color','black')
    hold on


    lgd = legend({'A=0.2(M=0.)','A=0.2(M=0.3)','A=0.2(M=0.5)',' A=0.5(M=0.)','A=0.5(M=0.3)', 'A=0.5(M=0.5)'},'FontSize',12,'TextColor','black','Location','NorthEast')

end


end
