function   my_nusselt
clc;close all;
c       =  {'-b',  '--b', '.:b'};
c1       =  {'-r',  '--r', '.:r'};
RC = .1:.3:1.5;
FW = .3:.3:2;
M = 1.5; BI = 1.5; BE = 1.5;R =.5;EC = 3; GR = 1.5; GAMMA = 5;  AE = 1 + BE * BI;
GM = 6; PR = .3  ;S = 0; S0 = 1; K = 1; SC = .6;
D = 2;EPS =.8; ST = .5; STT = .5;
NEBLA = 2; LAMB = .8; short = AE ^ 2 + BE ^ 2;
short=M/(AE^2+BE^2);
function ysol = bvpex1(x,y)        %bvpex1(ind variable, dependent variable)

yy3=-(D*y(7)+GR*y(8)+GM*y(10)+(y(1)*y(3))/(EPS^2)-(K+AE*short)*y(2)+(R-short*BE)*y(4)-GAMMA*y(2)^2)*EPS/(1+D);
yy5=-(y(1)*y(5)/EPS^2-(K+short*AE)*y(4)-(R-short*BE)*y(2)-GAMMA*y(4)^2)*EPS/(1+D);
yy7=-(y(2)*y(6)+y(1)*y(7)-2*LAMB*y(6)-LAMB*y(3))/NEBLA;
yy9=-((1+D)*PR*EC*(y(3)^2+y(5)^2)+short*EC*PR*(y(2)^2+y(4)^2)+PR*y(1)*y(9)-ST*PR*y(2));
yy11=-(SC*y(1)*y(11)-STT*SC*y(2)+S0*SC*yy9+RC(i)*SC*y(10));
ysol=[y(2);y(3);yy3;y(5);yy5;y(7);yy7;y(9);yy9;y(11);yy11];
end

% Here I define residual of boundary conditions
function res = bcex1(y0, yinf)


res=[y0(1)-FW(j);y0(2)-1;y0(4);y0(6)+y0(3)/2;y0(8)-1+ST/2;y0(10)-1+STT/2;
    yinf(2);yinf(4);yinf(6);yinf(8);yinf(10)];

end
fdp=zeros(6,1);
for i=1:length(RC)
    for j=1:length(FW)

        % Initial values
        sol1=bvpinit(linspace(0,12,600),[ 0 0 0 0 0 0 0 0 0 0 0]);
        % sol1 = bvpinit(linspace(xa,xb,n), [initiall guesses for f f' f'' etc]);

        % bvp4c takes two functions defined below and give solution in structure
        % form
        sol = bvp5c(@bvpex1, @bcex1, sol1);
        sol = bvp5c(@bvpex1, @bcex1, sol);
   

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
    plot(FW, fdp,c{i}, 'linewidth', 2.5);

    xlabel('Fw')
    ylabel('f''"(0)',Rotation=0)
    % title('Pr = 6.8, Rd=Kr=0.2, Nb=Nt=s=0.1, Le=Lb=Pe=0.5, Ec=0','fontsize',10,'fontweight','bold','color','black')
    
    hold on


end
title('Skin Fraction','FontWeight','bold',FontSize=15)
k=['Rc=',num2str(RC(1));'Rc=',num2str(RC(2));'Rc=',num2str(RC(3))];
lgd=legend(k);
hold off
end
% 2*f''(0)/Re_x         skin
% -theta'/sqrt(Re_x)       nusselt