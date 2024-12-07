clc;clear;
hold on
% sol=jobs(2,200);
% sol=jobs(2,50);
% sol=jobs(2,100);
sol=jobs(2,1000);
x=sol.x;
y=sol.y;
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
xlim([-.5 3]);
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

function out=jobs(ROW,LENGTH )
% RC=.1;FW=.8;
% M=0.5;AE = 2.1;BI = 0.7;BE = 0.4;R=0.6;GR=10.0;GM=5.0;PR=0.71;EC=0.01;S1 = 0.5;
% S0 = 1.0;K=0.5;SC=0.6;D = 0.5;BAR = 0.5;EPS = 0.6;ST = 0.5;STT= 0.5;LAMB = 3.0;
% short=AE^2+BE^2;
RC=.1;FW=.5;
M=0.4;BI = 0.2;BE = 0.1;AE = 1+BE*BI;R=0.6;GR=4.0;GM=2.0;PR=0.71;EC=0.1;S1 = 0.5;
S0 = .20;K=0.5;SC=0.6;D = 0.5;BAR = 0.5;EPS = 0.3;ST = 0.4;STT= 0.3;LAMB = 2.0;
short=AE^2+BE^2;XX=1+D;
    %options = bvpset('RelTol',0.0000001,'Stats','on')
    sol1=bvpinit(linspace(0,4,1000),[0 0 0 0 0 0 0 0 0 0 0 0 0]);
    sol=bvp4c(@bvp2d,@bc2d,sol1);
%     x=sol.x;
%     y=sol.y;
    out=sol;
 
    hold on
%     plot(x,y(2,:),'r');
%      plot(x,y(4,:),'g');
%     plot(x,y(6,:),'b');
%      plot(x,y(8,:),'c');
%     plot(x,y(10,:),'m');
%     plot(x,y(12,:),'y');
%     legend("2","4","6","8","10","12");
    
    %plot(x,sol.res)
    
    
    hold off
    
   % ylabel("Primary Velocity(F)")
    
    
    %xline(0);
    %yline(0);
    
    function yvector =bvp2d(t,y)
   
        yy3=-(-D*y(9)+GR*y(10)+GM*y(12)+y(1)*y(3)/(EPS^2)-...
            (K+M*AE/short)*y(2)+(R-(M*BE/short))*y(4)-BAR*y(2)^2)/((1+D)/EPS);
        yy5=-(D*y(7)+y(1)*y(5)/EPS^2-(K+M*AE/short)*y(4)-...
            (R-M*BE/short)*y(2)-BAR*y(4)^2)/((1+D)/EPS);
        yy7=-(y(6)*y(2)+y(1)*y(7))/LAMB;
        yy9=-(y(8)*y(2)+y(1)*y(9))/LAMB;
        Holder=-((1+D)*PR*EC*(y(3)^2+y(5)^2)+EC*M*PR*(y(2)^2+...
            y(4)^2)/short+PR*y(1)*y(11)-ST*PR*y(2));
        yy11=Holder;
        yy13=-(y(1)*SC*y(13)-STT*SC*y(2)+S0*SC*Holder-RC*SC*y(12));
        yvector=[y(2);y(3);yy3;y(5);yy5;y(7);yy7;y(9);yy9;y(11);yy11;y(13);yy13];
    end
    function residual=bc2d(y0,yinf)
    
        residual=[y0(2)-1;y0(1)-FW;y0(4);y0(6)+S1*y0(5);y0(8)-S1*y0(3);...
            y0(10)-1+ST/2;y0(12)-1+STT/2;...
            yinf(2);yinf(4);yinf(6);yinf(8);yinf(10);yinf(12)];
    end
%return sol;    
end