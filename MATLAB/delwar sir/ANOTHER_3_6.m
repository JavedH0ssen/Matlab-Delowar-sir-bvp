%% Final Shortened Complex one
clc; clear;
% sol=jobs(2,50);
% sol=jobs(2,100);
aa = .1;
bb = .01;
cc =.5;
name = 'Fw';
main(aa, bb, cc, name)

function main(aa, bb, cc, name)
    sol = jobs(aa);
    sol1 = jobs(bb);
    sol2 = jobs(cc);
        figure('Name',name)
    hold on
    subplot(2,2,1)
    plotting(2,"f'",sol,sol1,sol2,aa,bb,cc)%,sol1,sol2,aa,bb,
    subplot(2,2,2)
    plotting(5,"g",sol,sol1,sol2,aa,bb,cc)
    subplot(2,2,3)
    plotting(8,"\Theta",sol,sol1,sol2,aa,bb,cc)
    subplot(2,2,4)
    plotting(10,"\phi",sol,sol1,sol2,aa,bb,cc)

    
 
end

function plotting(line, titles, sol, sol1, sol2, aa, bb, cc, name)
  
    hold on
    a = plot(sol.x, sol.y(line, :), '-r');
    b = plot(sol1.x, sol1.y(line, :), '--g');
    c = plot(sol2.x, sol2.y(line, :), '-.b');
    xline(0);
    yline(0);
%     lgd = legend([a; b; c], {append(name, '=', num2str(aa)); ...
%         append(name, '=', num2str(bb)); append(name, '=', num2str(cc))}, 'Location', 'northeast');
%         lgd.FontSize = 7;  % Set the font size to 14
    xline(0);
    yline(0);
 
    box on
    hold off

end

function out=jobs(v)
RC=.5;FW=.6;
M=0.5;BI =5;BE = 1.5;AE = 1+BE*BI;R=0.6;GR=8;
GM=5;PR=.01;EC=v; S = 0.5;S0 =.1;K=.5;SC=.6;
D = 2;GAMMA = 1;EPS = .8;ST = 5;STT= .5;
NEBLA = 2;LAMB=.8;short=AE^2+BE^2;XX=1+D;W1=1.5;CP=1.5;DF=.8;n=.01;
  disp('1');
  aaa=3.5;
  sol1 = bvpinit(linspace(0,aaa,500),[1 0 1 1 0 1 1 1 1 1 1]);
sol = bvp4c(@bvp2d, @bc2d, sol1); 
  disp('2');
% sol1=bvpinit(sol,[0 1.1]);
% sol = bvp4c(@bvp2d, @bc2d, sol1);
sol1=bvpinit(sol,[0 aaa]);
sol = bvp4c(@bvp2d, @bc2d, sol1);
sol1=bvpinit(sol,[0 aaa]);
sol = bvp4c(@bvp2d, @bc2d, sol1);
% %   disp('3');
% sol1=bvpinit(sol,[1.5 1.87]);
% sol = bvp4c(@bvp2d, @bc2d, sol1);
%   disp('4');
% sol1=bvpinit(sol,[0 1.85]);
% sol = bvp4c(@bvp2d, @bc2d, sol1);
%   disp('5');
% sol1=bvpinit(sol,[0 1.86]);
% sol = bvp4c(@bvp2d, @bc2d, sol1);
%   disp('6');
% sol1=bvpinit(sol,[1.87 1.89]);
% sol = bvp4c(@bvp2d, @bc2d, sol1);
%   disp('7');
% sol1=bvpinit(sol,[1.89 1.9]);
% sol = bvp4c(@bvp2d, @bc2d, sol1);
%   disp('8');
% sol1=bvpinit(sol,[1.9 2.1]);
% sol = bvp4c(@bvp2d, @bc2d, sol1);
%   disp('9');
% sol1=bvpinit(sol,[2 2.5]);
% sol = bvp4c(@bvp2d, @bc2d, sol1);
%   disp('10');
% sol1=bvpinit(sol,[2.5 3]);
% sol = bvp4c(@bvp2d, @bc2d, sol1);
%   disp('11');
% sol1=bvpinit(sol,[0 3]);
% sol = bvp4c(@bvp2d, @bc2d, sol1);
%   disp('12');
% sol1=bvpinit(sol,[0 11.3]);
% sol = bvp4c(@bvp2d, @bc2d, sol1);
%   disp('13');
% sol1=bvpinit(sol,[0 15.5]);
% sol = bvp4c(@bvp2d, @bc2d, sol1);
    fprintf('################################################')
%     fprintf(sol.stats.maxres)
    a=sol.stats.maxres
    fprintf('################################################')
    out=sol;
    function yvector =bvp2d(~,y)
%         if y(1)==0
%             y(1)=0.01;
%         end
        yy4=(-y(4)-y(3)*y(1)-W1*(2*y(4)*y(2)-y(3)^2)+K*y(2)+(M/(AE^2+BE^2))*...
            (AE*y(2)+BE*y(5))-R*y(5)+GAMMA*y(2)^2-GM*y(10)-GR*y(8))/(W1*y(1));
        yy7=(-y(1)*y(6)-y(7)-2*W1*y(2)*y(7)+W1*y(6)^2+K*y(5)+(M/short)*(AE*y(5)-...
            BE*y(2))+R*y(2)+GAMMA*y(5)^2)/(W1*y(1));
%         yy11=(y(1)*y(9)-y(1)*y(11)*DF*SC+(EC/CP)*(y(3)^2+y(6)^2)+...
%             (M*EC)/short*(y(2)^2+y(5)^2)+RC-DF*SC*(y(10)+1))/(1/PR-DF*S0*SC);

%%         %one taking Df=0; \phi_2 in eqn 3 result in y3 and y6,7
        yy9=(y(1)*y(9)+EC*(y(3)^2+y(6)^2)+M*EC/short*(y(2)^2+y(5)^2))*(-PR);
        yy11=(-y(1)*y(11)-S0*yy9+RC*y(10))*SC;
       

%%       %one taking S0 zero \theta in eqn 4
%         yy11=(-y(1)*y(11)+RC*y(10))*SC;
%         yy9=(y(1)*y(9)+EC*(y(3)^2+y(6)^2)+M*EC/short*(y(2)^2+y(5)^2)+DF*yy11)*(-PR);
 
%%        one torsten suggested
%         yk9=(y(1)*y(9)+EC*(y(3)^2+y(6)^2)+M*EC/short*(y(2)^2+y(5)^2))*PR;
%         yk11=(-y(1)*y(11)+RC*y(10));  
%         yy11=(yk9*S0+yk11)/-(PR*DF*S0-1/SC);
%         yy9=(yk9+yk11*PR*DF*SC)/-(1-SC*S0*PR*DF);
yvector=[y(2);y(3);y(4);yy4;y(6);y(7);yy7;y(9);yy9;y(11);yy11];
    end
    function residual=bc2d(y0,yinf)
    
        residual=[y0(1)-FW;y0(2)-1;y0(5);y0(8)-1;y0(10)-1;yinf(3);yinf(6);
            yinf(2);yinf(5);yinf(8);yinf(10)];
    end
   
end
