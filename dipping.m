clc
close all;clear;
tic;
current="\Lambda";
a=1;
b=1.5;
c=.7;
sol1=jobs(a);
% sol2=jobs(b);
% sol3=jobs(c);
sol2=sol1;
sol3=sol1;
size=500;


main(sol1,sol2,sol3,size,a,b,c,current);
toc
function main(sol1,sol2,sol3,size,a,b,c,current)
figure;
hold on
n=1;
% f=figure(n);
% f.Position = [200 200 size size];
subplot(2,2,1);
ploting("f'",sol1,sol2,sol3,2,current,a,b,c,n);
% xlim([0 8])
n=n+1;
% f=figure(n);
subplot(2,2,2)
% f.Position = [100 100 size size];
ploting("g",sol1,sol2,sol3,5,current,a,b,c,n);
legend('Location', 'southeast');
n=n+1;
% f=figure(n);
subplot(2,2,3)
% f.Position = [100 100 size size];
ploting("\theta",sol1,sol2,sol3,8,current,a,b,c,n);
n=n+1;
% f=figure(n);
% f.Position = [100 100 size size];
subplot(2,2,4)
ploting("\phi",sol1,sol2,sol3,10,current,a,b,c,n);
n=n+1;
% f=figure(n);
% f.Position = [100 100 size size];
% ploting("f",sol1,sol2,sol3,1,current,a,b,c,n)
% n=n+1;

hold off
end

function a=ploting(title,sol1,sol2,sol3,nn,cur,a,b,c,n)
hold on;
plot(sol1.x,sol1.y(nn,:),'r');
plot(sol2.x,sol2.y(nn,:),'g');
plot(sol3.x,sol3.y(nn,:),'b');
ylabel(title,'FontWeight', 'bold',Rotation=0);
xlabel("\eta",'FontWeight', 'bold');
x=.028;
y=.028;
annotation("arrow", [0.1936 0.5], [x x]);
annotation("arrow", [y y], [0.1646 0.48]);
set(gca,"FontSize",15);
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
RC = 2; FW = 1.5;
M = 1.5; BI = 1.5; BE = 1.5;R =.5;EC = 15; GR = 1.5; GAMMA = .5;  AE = 1 + BE * BI;
GM = .6; PR = .3  ;S = 0; S0 = .5; K = 1; SC = .6;
D = 2;EPS =.8; ST = .5; STT = .5;Wi=.5;
NEBLA = 2; LAMB = .8; short = AE ^ 2 + BE ^ 2;DF=.5;
short=M/(AE^2+BE^2);
%options = bvpset('RelTol',0.0000001,'Stats','on')
sol1=bvpinit(linspace(0,10,3000),[ FW 1 0 0 1 0 1 0 0 1 1]);
sol1=bvp4c(@bvp2d,@bc2d,sol1);
n=0;
while sol1.stats.maxres>1e-6 && n<15
    sol1=bvp4c(@bvp2d,@bc2d,sol1);
    n=n+1;
    % if n>10 && sol1.stats.maxres>2
    %     break;
    % end
end
fprintf(' %d     %d \n',n,sol1.stats.maxres)

sol=bvp4c(@bvp2d,@bc2d,sol1);

out=sol;




    function yvector =bvp2d(~,y)
        yy4=(y(4)+y(1)*y(3)-Wi*(-2*y(2)*y(4)+(y(3)^2))-K*y(2)-short*(AE*y(2)+BE*y(5))+GR*y(8)+GM*y(10)+R*y(5)-GAMMA*y(2)^2)/(-y(1)*Wi);
        yy7=(y(7)+y(1)*y(6)-Wi*(-2*y(2)*y(7)+y(6)^2)-K*y(5)+short*(BE*y(2)-AE*y(5))-R*y(2)-GAMMA*y(5)^2)/(-Wi*y(1));
        
        theta=-(PR*y(1)*y(9)+PR*EC*(y(3)^2+y(6)^2)+PR*EC*short*(y(2)^2+y(5)^2));  %-pr*df*phi'' 
        phi=-(SC*y(1)*y(11)-RC*SC*y(10));% -SC*S0*theta

%%        taking df=0 evaluate theta
        yy9=theta;
        yy11=phi-SC*S0*yy9;
%%      taking S0=0 evaluate phi
        % yy11=phi;
        % yy9=theta-PR*DF*yy11;
%%      taking both yy9
        % yy9=(theta+PR*DF*phi)/(1+SC*S0*PR*DF);
        % yy11=phi-SC*S0*yy9;
%%    taking both retry 
        % yy11=(phi+theta*SC*S0)/(1+PR*DF*SC*S0);
        % yy9=theta-PR*DF*yy11;
        
        yvector=[y(2);y(3);y(4);yy4;y(6);y(7);yy7;y(9);yy9;y(11);yy11;];

    end
    function residual=bc2d(y0,yinf)


        residual=[ y0(1)-FW;y0(2)-1; y0(5); y0(8)-1; y0(10)-1; y0(3);yinf(3);
           yinf(2); yinf(5); yinf(8); yinf(10) ];
    end
%return sol;
end