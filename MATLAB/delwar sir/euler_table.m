clc; clearvars; close all;
f = @(x, y) (y - x) / (x + y);
x0 = 0;
xn = 0.1;
h = 0.025;
x = x0:h:xn;
y0 = 1;
y = zeros(size(x));
y(1) = y0;
for i = 1:length(x) - 1
y(i+1) = y(i) + h*f(x(i), y(i));
end
T = table(x', y', 'VariableNames', {'x', 'y'});
disp(T);