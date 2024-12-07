%% Final Shortened Complex one
clc; clear;
% warning('off','all')
aa = 1.5;
bb = 3;
cc =2;
% aa = .01;
% bb = .1;
% cc =1;
name = 'Rc';
tic;
main(aa, bb, cc, name)
toc
function main(aa, bb, cc, name)
    sol = jobs(aa);
    sol1 = jobs(bb);
    sol2 = jobs(cc);
        figure('Name',name)
    hold on
    subplot(2,2,1)
    plotting(2,"f'",sol,sol1,sol2,aa,bb,cc,name)%,sol1,sol2,aa,bb,
    subplot(2,2,2)
    plotting(5,"g",sol,sol1,sol2,aa,bb,cc,name)
    subplot(2,2,3)
    plotting(8,"\Theta",sol,sol1,sol2,aa,bb,cc,name)
%     xlim([0 4.5])
    subplot(2,2,4)
    plotting(10,"\phi",sol,sol1,sol2,aa,bb,cc,name)
%     xlim([0 3.5])
    
 
end

function plotting(line, titles, sol, sol1, sol2, aa, bb, cc, name)
  
    hold on
    xline(0);yline(0);
    a = plot(sol.x, sol.y(line, :), '-r');
    b = plot(sol1.x, sol1.y(line, :), '--g');
    c = plot(sol2.x, sol2.y(line, :), '-.b');
 lgd=legend([a; b; c], {[name '= ' num2str(aa)],[name '= ' num2str(bb)],[name '= ' num2str(cc)]});

%     lgd = legend([a; b; c], {append(name, '=', num2str(aa)); ...
%         append(name, '=', num2str(bb)); append(name, '=', num2str(cc))}, 'Location', 'northeast');
%         lgd.FontSize = 7;  % Set the font size to 14
%     xline(0);
%     yline(0);
    xlabel('\eta','FontWeight','bold','FontSize', 13, 'Interpreter', 'tex')
    ylabel([titles],'FontWeight','bold','FontSize', 13, 'Interpreter', 'tex')
    hYLabel = get(gca,'YLabel');
    set(hYLabel,'Units','normalized');
    set(hYLabel,'rotation',0)
    set(hYLabel,'Position',get(hYLabel,'Position') + [-0.02 0 0]);
%     xlim([0 5])
    box on
    
    hold off

end

function out=jobs(v)
RC=1.5;FW=v;%fw=1.89
M=7;BI =1.5;BE = 1.5;AE = 1+BE*BI;R=20;GR=1.35;
GM=3;PR=.71;EC=.3;S0 = .1;K=.1;SC=.6;
GAMMA = 1;short=AE^2+BE^2;W1=3.5;DF=.5;
AAA=8;
  sol1 = bvpinit(linspace(0,AAA,500),[1 0 1 1 0 1 1 1 1 1 1]);
sol = bvp4c(@bvp2d, @bc2d, sol1); 
% sol1=bvpinit(sol,[0 1.1]);
% sol = bvp4c(@bvp2d, @bc2d, sol1);
% sol1=bvpinit(sol,[0 .5]);
% sol = bvp4c(@bvp2d, @bc2d, sol1);
n=0;
  while(sol.stats.maxres>1e-5)
      lo=sol.stats.maxres;
      sol1=bvpinit(sol,[0 AAA]);
      sol = bvp4c(@bvp2d, @bc2d, sol1);
      if ((sol.stats.maxres==lo && n>2)||n>50)
          break;
      end
      n=n+1;
  end
fprintf('done (%d) = %e\n',n,sol.stats.maxres);
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
%         yy9=(y(1)*y(9)+EC*(y(3)^2+y(6)^2)+M*EC/short*(y(2)^2+y(5)^2))*(-PR);
%         yy11=(-y(1)*y(11)-S0*yy9+RC*y(10))*SC;
       

%%       %one taking S0 zero \theta in eqn 4
%         yy11=(-y(1)*y(11)+RC*y(10))*SC;
%         yy9=(y(1)*y(9)+EC*(y(3)^2+y(6)^2)+M*EC/short*(y(2)^2+y(5)^2)+DF*yy11)*(-PR);
 
%%        one torsten suggested
        yk9=(y(1)*y(9)+EC*(y(3)^2+y(6)^2)+M*EC/short*(y(2)^2+y(5)^2))*PR;
        yk11=(-y(1)*y(11)+RC*y(10));  
        yy11=(yk9*S0+yk11)/-(PR*DF*S0-1/SC);
        yy9=(yk9+yk11*PR*DF*SC)/-(1-SC*S0*PR*DF);
yvector=[y(2);y(3);y(4);yy4;y(6);y(7);yy7;y(9);yy9;y(11);yy11];
    end
    function residual=bc2d(y0,yinf)
    
        residual=[y0(1)-FW;y0(2)-1;y0(5);y0(8)-1;y0(10)-1;yinf(3);yinf(6);
            yinf(2);yinf(5);yinf(8);yinf(10)];
    end
   
end

