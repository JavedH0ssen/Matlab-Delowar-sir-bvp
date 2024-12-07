k=jobs(1);
eta = linspace (0, 6, 100);
figure
subplot(2,2,1)
plot(eta,k(3,:));
subplot(2,2,2)
plot(eta,k(5,:));
function out=jobs(NDEX)
% RC=.1;FW=.8;
% M=0.5;AE = 2.1;BI = 0.7;BE = 0.4;R=0.6;GR=10.0;GM=5.0;PR=0.71;EC=0.01;S1 = 0.5;
% S0 = 1.0;K=0.5;SC=0.6;D = 0.5;BAR = 0.5;EPS = 0.6;ST = 0.5;STT= 0.5;LAMB = 3.0;
% short=AE^2+BE^2;
RC=.1;FW=0.5;
M=0.5;BI = 0.3;BE = 0.4;AE = 1+BE*BI;R=0.6;GR=4;GM=2;PR=0.71;EC=NDEX; S = 0.5;
S0 = 1;K=0.5;SC=0.6;D = 2;BAR = 0.5;EPS = 0.6;ST = 0.5;STT= 0.5;NEBLA = 2.0;LAMB=.5;
short=AE^2+BE^2;XX=1+D;
    %options = bvpset('RelTol',0.0000001,'Stats','on')
    sol1=bvpinit(linspace(0,6,100),[0 0 0 0 0 0 0 0 0 0 0]);
    sol=bvp4c(@bvp2d,@bc2d,sol1);
%     x=sol.x;
%     y=sol.y;
    eta = linspace (0, 6, 100);
    y= deval (sol,eta);
    fprintf('f"(0)=%7.9f(reduced skin friction)\n',y(3));  % reduced skin friction 
    fprintf('-theta(0)=%7.9f(reduced local nusselt number)\n',-y(5)); %reduced local nusselt number
    fprintf('Cfx(Skin Fraction)=%7.9f\n',(y(3)));    %  skin friction 
    fprintf('Nux(Nusselt Number=%7.9f\n)',-y(5));     % local nusselt number
    fprintf('#############################\n');
    out = y;
    function yvector =bvp2d(~,y)
   
        yy3=-(D*y(7)+GR*y(8)+GM*y(10)+y(1)*y(3)/(EPS^2)-...
            (K+M*AE/short)*y(2)+(R-(M*BE/short))*y(4)-BAR*y(2)^2)/((1+D)/EPS);
        yy5=-(y(1)*y(5)/EPS^2-(K+M*AE/short)*y(4)-...
            (R-M*BE/short)*y(2)-BAR*y(4)^2)/((1+D)/EPS);
        yy7=-(y(6)*y(2)+y(1)*y(7)-2*LAMB*y(6)-LAMB*y(3))/NEBLA;
        Holder=-((1+D)*PR*EC*(y(3)^2+y(5)^2)+EC*M*PR*(y(2)^2+...
            y(4)^2)/short+PR*y(1)*y(9)-ST*PR*y(2));
        yy9=Holder;
        yy11=-(y(1)*SC*y(11)-STT*SC*y(2)+S0*SC*Holder+RC*SC*y(10));
        yvector=[y(2);y(3);yy3;y(5);yy5;y(7);yy7;y(9);yy9;y(11);yy11];
    end
    function residual=bc2d(y0,yinf)
    
        residual=[y0(2)-1;y0(1)-FW;y0(4);y0(6)+y0(3)/2;y0(8)-1+ST/2;...
            y0(10)-1+STT/2;...
            yinf(2);yinf(4);yinf(6);yinf(8);yinf(10)];
    end
%return sol;    
end