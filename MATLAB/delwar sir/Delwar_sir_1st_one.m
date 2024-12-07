%% First Trial one from Delwar Sir
clear;clc;
%PR= 1.38;
%AE = 1.4; BE = 0.5;M= 1.2;G =0.2;EC=0.1; R = 0.6; FW = 1.5;GR = 10.0;
figure

ROW=2;
subplot(2,2,1)
sbplot(2,3)
title("Primary Velocity(F) ")
subplot(2,2,2)
sbplot(4,3.5)
title("Secondary Velocity(G) ")
subplot(2,2,3)
sbplot(6,2.5)
title("Temparature (\theta) ")
function sbplot(ROW,LENGTH)
    hold on
    jobs(1.38,ROW,LENGTH);
    jobs(1,ROW,LENGTH);
    jobs(0.7,ROW,LENGTH);
    legend("1.38","1","0.7")
    hold off;
end
function jobs(PR,ROW,LENGTH )

AE = 1.4; BE = 0.5;M= 1.2;G =0.2;EC=0.1; R = 0.6; FW = 1.5;GR = 10.0;
    sol1=bvpinit(linspace(0,LENGTH,100),[0 0 0 0 0 0 0]);
    sol=bvp4c(@bvp2d,@bc2d,sol1);
    x=sol.x;
    y=sol.y;
    %legend('Plot',string(PR))
    %hold on
    plot(x,y(ROW,:));%,x,y(4,:),"g";x,y(6,:),"b"2
    
    
    %hold off
    
   % ylabel("Primary Velocity(F)")
    
    
    %xline(0);
    %yline(0);
    %legend();
    function yvector =bvp2d(~,y)
        %PR= 1.38;
   % AE = 1.4; BE = 0.5;M= 1.2;G =0.2;EC=0.1; R = 0.6; FW = 1.5;GR = 10.0;
        X=(G+M*AE/(AE^2+BE^2));
        Z=(R-M*BE/(AE^2+BE^2));
        yy1=-(y(1)*y(3)-X*y(2)+Z*y(4)+GR*y(6));
        yy2=-(y(1)*y(5)-X*y(4)-Z*y(2));
        yy3=-(PR*y(1)*y(7)-2*PR*y(2)*y(6)+PR*EC*(y(3)^2+y(5)^2)...
            +PR*EC*M*(y(2)^2+y(4)^2)/(AE^2+BE^2));
    
        yvector=[y(2);y(3);yy1;y(5);yy2;y(7);yy3];
    end
    function residual=bc2d(y0,yinf)
    %FW=1.5;
        residual=[y0(2)-1;y0(1)-FW;y0(4);y0(6)-1;yinf(2);yinf(4);yinf(6)];
    end
end