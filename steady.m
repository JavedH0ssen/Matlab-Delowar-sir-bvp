clc;clear;
% for n=1:5
%     close(figure(n));
% end
close all;
current="\Lambda";
a=.1;
b=.5;
c=.8;
sol1=jobs(a);
sol2=jobs(b);
sol3=jobs(c);
size=500;


main(sol1,sol2,sol3,size,a,b,c,current)
function main(sol1,sol2,sol3,size,a,b,c,current)
figure
hold on
n=1;
f=figure(n);
f.Position = [200 200 size size];
ploting("f",sol1,sol2,sol3,2,current,a,b,c,n)
% xlim([0 8])
n=n+1;
f=figure(n);
f.Position = [100 100 size size];
ploting("g",sol1,sol2,sol3,3,current,a,b,c,n)
legend('Location', 'southeast');
n=n+1;
f=figure(n);
f.Position = [100 100 size size];
ploting("\Xi",sol1,sol2,sol3,6,current,a,b,c,n)
n=n+1;
f=figure(n);
f.Position = [100 100 size size];
ploting("\theta",sol1,sol2,sol3,8,current,a,b,c,n)
n=n+1;
f=figure(n);
f.Position = [100 100 size size];
ploting("\phi",sol1,sol2,sol3,10,current,a,b,c,n)
n=n+1;

hold off
end

function a=ploting(title,sol1,sol2,sol3,nn,cur,a,b,c,n)
hold on
plot(sol1.x,sol1.y(nn,:),'r');
plot(sol2.x,sol2.y(nn,:),'g');
plot(sol3.x,sol3.y(nn,:),'b');
ylabel(title,'FontWeight', 'bold',Rotation=0)
xlabel("\eta",'FontWeight', 'bold');
x=.028;
y=.028;
annotation("arrow", [0.1936 0.5], [x x]);
annotation("arrow", [y y], [0.1646 0.48]);
set(gca,"FontSize",15)
legend([cur + " = " + num2str(a);cur + " = " ...
    + num2str(b); cur + " = " + num2str(c)]);
% set(gcf,'color','white')
box on
yadd=[sol1.y sol2.y sol3.y];
ymax= max(yadd(nn,:));
ymin=min(yadd(nn,:));
ylim([ymin ymax]);

xlim_vals = xlim;
ylim_vals = ylim;
% set(gca, 'XTick', linspace(xlim_vals(1), xlim_vals(2), 5));
% set(gca, 'YTick', round(linspace(ylim_vals(1), ylim_vals(2), 10), 2));
% pbaspect([1 1 1])
tightInset = get(gca, 'TightInset');
extraSpace = 0.01;
looseInset = tightInset + extraSpace;
set(gca, 'LooseInset', looseInset);
name_fig=strcat(cur+num2str(n)+".png");
% saveas(gcf,name_fig)
% ax = gca;
% exportgraphics(ax,name_fig,'Resolution',500)

hold off
end

function out=jobs(v)
RC = v; FW = 1.5;
M = 1.5; BI = 1.5; BE = 1.5;R =.5;EC = 3; GR = 1.5; GAMMA = 5;  AE = 1 + BE * BI;
GM = 6; PR = .3  ;S = 0; S0 = 1; K = 1; SC = .6;
D = 2;EPS =.8; ST = .5; STT = .5;
NEBLA = 2; LAMB = .8; short = AE ^ 2 + BE ^ 2;
short=M/(AE^2+BE^2);
%options = bvpset('RelTol',0.0000001,'Stats','on')
sol1=bvpinit(linspace(0,11,1000),[ 0 0 .1 0 0 0 0 0 0 0 0]);


sol1=bvp4c(@bvp2d,@bc2d,sol1);
sol1=bvp4c(@bvp2d,@bc2d,sol1);
sol=bvp4c(@bvp2d,@bc2d,sol1);
out=sol;


    function yvector =bvp2d(~,y)
        yy3=-(D*y(7)+GR*y(8)+GM*y(10)+(y(1)*y(3))/(EPS^2)-(K+AE*short)*y(2)+(R-short*BE)*y(4)-GAMMA*y(2)^2)*EPS/(1+D);
        yy5=-(y(1)*y(5)/EPS^2-(K+short*AE)*y(4)-(R-short*BE)*y(2)-GAMMA*y(4)^2)*EPS/(1+D);
        yy7=-(y(2)*y(6)+y(1)*y(7)-2*LAMB*y(6)-LAMB*y(3))/NEBLA;
        yy9=-((1+D)*PR*EC*(y(3)^2+y(5)^2)+short*EC*PR*(y(2)^2+y(4)^2)+PR*y(1)*y(9)-ST*PR*y(2));
        yy11=-(SC*y(1)*y(11)-STT*SC*y(2)+S0*SC*yy9-RC*SC*y(10));
        
        yvector=[y(2);y(3);yy3;y(5);yy5;y(7);yy7;y(9);yy9;y(11);yy11];
    end
    function residual=bc2d(y0,yinf)


        residual=[y0(1)-FW;y0(2)-1;y0(4);y0(6)+y0(3)/2;y0(8)-1+ST/2;y0(10)+1+STT/2;
            yinf(2);yinf(4);yinf(6);yinf(8);yinf(10)];
    end
%return sol;
end