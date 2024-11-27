aa=jobs(1)
function out=jobs(v)
RC = .1; FW = 1.5;
M = 1.5; BI = 1.5; BE = 1.5;R =.5;EC = 3; GR = 1.5; GAMMA = 5;  AE = 1 + BE * BI;
GM = 6; PR = .3  ;S = 0; S0 = 1; K = 1; SC = .6;
D = 2;EPS =.8; ST = .5; STT = .5;Wi=.5;
NEBLA = 2; LAMB = .8; short = AE ^ 2 + BE ^ 2;DF=.5;
short=M/(AE^2+BE^2);
%options = bvpset('RelTol',0.0000001,'Stats','on')
sol1=bvpinit(linspace(0,6,1000),[ 0 0 0 0 0 0 0 0 0 0 0]);
sol=bvp4c(@bvp2d,@bc2d,sol1);
out=sol;




    function yvector =bvp2d(~,y)
        yy4=(y(4)+y(1)*y(3)-Wi*(-2*y(2)*y(4)+(y(3)^2))-k*y(2)-short*(AE*y(2)+BE*y(5))+GR*y(8)+GM*y(10)+R*y(5)-GAMMA*y(2)^2)/(-y(1)*Wi);
        yy7=(y(7)+y(1)*y(6)-Wi*(-2*y(2)*y(7)+y(6)^2)-k*y(5)+short*(BE*y(2)-AE*y(5))-R*y(2)-GAMMA*y(5)^2)/(-Wi*y(1));
        
        theta=-(PR*y(1)*y(9)+PR*EC*(y(3)^2+y(6)^2)+PR*EC*short*(y(2)^2+y(5)^2));  %-pr*df*phi'' 
        phi=-(SC*y(1)*y(11)-RC*SC*y(10));% -SC*S0*theta

%%        taking df=0 evaluate theta
        yy9=theta;
        yy11=phi-SC*S0*yy9;
%%      taking S0=0 evaluate phi
%         yy11=phi;
%         yy9=theta-PR*DF*yy11;
% %%      taking both yy9
%         yy9=(theta+PR*DF*phi)/(1+SC*S0*PR*DF);
%         yy11=phi-SC*S0*yy9;
% %%    taking both retry 
%         yy11=(phi+theta*SC*S0)/(1+PR*DF*SC*S0);
%         yy9=theta-PR*DF*yy11;
        
        yvector=[y(2);y(3);y(4);yy4;y(6);y(7);yy7;y(9);yy9;y(11);yy11;];

    end
    function residual=bc2d(y0,yinf)


        residual=[y0(2)-1;y0(1)-FW;y0(5);y0(8)-1;y0(10)-1;
           yinf(2);yinf(5);yinf(8);yinf(10) ];
    end
%return sol;
end