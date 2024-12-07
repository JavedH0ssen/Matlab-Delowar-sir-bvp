%% 3 Variable Failed One
clc;clear;
% sol=jobs(2,200);
% sol=jobs(2,50);
% sol=jobs(2,100);
sol=jobs(2,500);
x=sol.x;
y=sol.y;
figure

hold on

subplot(3,2,1)
plot(x,y(2,:));
title("Pimary Velocity(F')")
xline(0);
yline(0);
subplot(3,2,2)
plot(x,y(4,:));
title("Secondary Velocity(G)")
xline(0);
yline(0);
subplot(3,2,3)
plot(x,y(6,:));
title("Concentation 1 (\Xi _{1})")
xline(0);
yline(0);
subplot(3,2,4)
plot(x,y(8,:));
title("Concentation 2 (\Xi _{2})")
xline(0);
yline(0);
subplot(3,2,5)
plot(x,y(10,:));
title("Temparature(\Theta)")
xline(0);
yline(0);
subplot(3,2,6)
plot(x,y(12,:));
title("Non Partical Volume Fraction(\Phi)")
xline(0);
yline(0);

    
hold off

function sol2=jobs(ROW,LENGTH )
RC=.1;FW=.5;
M=0.4;BI = 0.2;BE = 0.1;AE = 1+BE*BI;R=0.6;GR=4.0;GM=2.0;PR=0.71;EC=0.1;S1 = 0.5;
S0 = .20;K=0.5;SC=0.6;D = 0.5;BAR = 0.5;EPS = 0.3;ST = 0.4;STT= 0.3;LAMB = 2.0;
short=AE^2+BE^2;XX=1+D;

    %options = bvpset('RelTol',0.0000001,'Stats','on')
%     options = bvpset('RelTol', 10e-4,'AbsTol',10e-7);
    sol1=bvpinit(linspace(0,5,500),[1 0 0 0 0 0 0 0 0 0 0 0 0]);
%     sol1=bvpinit(sol,[10 20]);                                  


    sol2=bvp4c(@bvp2d,@bc2d,sol1);
    sol1=bvpinit(sol2,[0 6]);
    sol3=bvp4c(@bvp2d,@bc2d,sol1);
%     sol1=bvpinit(sol3,[0 15]);
%     sol3=bvp4c(@bvp2d,@bc2d,sol1);
%     sol1=bvpinit(sol3,[0 20]);
%     sol3=bvp4c(@bvp2d,@bc2d,sol1);
    
    sol2=sol3;
%     
   
%     a(2)=plot(x,y(2,:));
%       a(4)=plot(x,y(4,:));
%     a(6)=plot(x,y(6,:));
%     a(8)= plot(x,y(8,:));
%    a(10)= plot(x,y(10,:));
%     plot(x,y(12,:));
%     legend("2","4","6","8","10","12");
    
    %plot(x,sol.res)
    
    
    
    
   % ylabel("Primary Velocity(F)")
    
    
    %xline(0);
    %yline(0);
    
    function yvector =bvp2d(t,y)
        yy3=EPS*D*y(9)/XX+EPS*GR*y(10)/XX-EPS*GM*y(12)/XX-(1/XX)*(1/EPS)*y(1)*y(3)+...
            (EPS/XX)*(K+M*AE/short)*y(2)-(EPS/XX)*(R-M*BE/short)*y(4);
        yy5=-(D*EPS/XX)*y(7)-(1/XX)*(1/EPS)*y(1)*y(5)+(EPS/XX)*(K+M*AE/short)*y(4)+...
            (EPS/XX)*(R-M*BE/short)*y(2)+(EPS/XX)*BAR*(y(2))^2;
        yy7=(-1/LAMB)*(y(6)*y(2)+y(1)*y(7));
        yy9=(-1/LAMB)+y(8)*y(2)+(-1/LAMB)*y(1)*(9);
        yy11=-XX*PR*EC*(y(3)^2+y(5)^2)-PR*y(1)*y(11)-(EC*M*PR/short)*(y(2)^2+y(4)^2)+ST*PR*y(2);
        yy13=-y(1)*SC*y(13)+STT*SC*y(2)-S0*SC*yy11+RC*SC*y(12);
        yvector=[y(2);y(3);yy3;y(5);yy5;y(7);yy7;y(9);yy9;y(11);yy11;y(13);yy13];
    end
    function residual=bc2d(y0,yinf)
    
        
    residual=[y0(2)-1;y0(1)-FW;y0(4);y0(6)+S1*y0(5);y0(8)-S1*y0(3);y0(10)-1+ST/2;y0(12)-1+STT/2;...
        yinf(2);yinf(4);yinf(6);yinf(8);yinf(10);yinf(12);];
    end
%return sol;    
end