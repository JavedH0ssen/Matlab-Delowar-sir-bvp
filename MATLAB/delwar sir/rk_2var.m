clc; clearvars; close all;
f = @(h, p) p*h + 1;
h0 = 0.0;
p0 = 1;
g = 0.2;
xn = 0.4;
n = (xn - h0) / g;
h(1) = h0; p(1) = p0;
for i = 1:n
h(i + 1) = h0 + i * g;
k1 = g * f(h(i), p(i));
k2 = g * f(h(i) + (g / 2), p(i) + (k1 / 2));
k3 = g * f(h(i) + (g / 2), p(i) + (k2 / 2));
k4 = g * f(h(i) + g, p(i) + k3);
p(i + 1) = p(i) + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4);
fprintf('Pressure at %0.1fm = %0.4f atm\n', h(i + 1), p(i + 1))
end