
function  Md_Yousuf_2024_Paper
clc
clear all
E = 0.003; Nt = 0.3; Rp = 1; Nb = 0.3; Le = 1; Kr = 0.20; AE = 2; M = 1; p = 1; Gr= 0.1; Gm = 0.1; theta_w = 0.5; Pr = 1.7;

lambda = 0.2; n = 0.3;

solinit  = bvpinit(linspace(0,10, 20), [1 0 0 0 0 0 0]);

sol      = bvp4c(@odefun, @bcfun, solinit );

x   = sol.x;  % domain

y   = sol.y  % y, y^' etc

skin_friction_coefficient = y(5,1)

figure (1)
plot(x, y(2, :))


    function dydx = odefun(x, y)

        temp1 = y(2)^2 - y(1)*y(3) + M*(y(2)*(sin(pi*x))^2 - E) - Gr*y(4) - Gm*y(6);

        temp2 = 1/(1 + Rp*(1 + (theta_w - 1)*y(4))^3)*(-Pr*y(5)*y(1) - Nb*y(5)*y(7) - Nt*y(5)^2 -3*Rp*(theta_w - 1)*(1 +(theta_w - 1)*y(4))^2*y(5)^2);

        temp3 = -Le*y(7)*y(1) - Nt/Nb*temp1 + Kr*Le*(1 + lambda*y(4))^n*exp(-AE/(1 + lambda*y(4)));

              dydx = [y(2); y(3); temp1; y(5); temp2; y(7); temp3];
    end




    function residual = bcfun(y0, yinf)

             residual = [y0(2) - 1; y0(1); y0(6) - 1; y0(4) - 1; yinf(2); yinf(4); yinf(6)];

    end
end