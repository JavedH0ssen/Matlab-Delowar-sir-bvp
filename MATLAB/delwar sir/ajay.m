close all
clc
clear
figure
hold on
subplot(2,2,1)
x=[0,6.028,6.868,7.318,8.108,14.154]
y=[25,24,23,22,20.5,22.5]
plot(x,y)
xlabel('Distance of section from section a (cm)');
ylabel('Piezometer head (cm)');
title('Hydraulic gradient line for observation 2')
grid on
subplot(2,2,2)
x=[0,6.028,6.868,7.318,8.108,14.154]
y=[22.5,22,21,20,18,20]
plot(x,y)
xlabel('Distance of section from section a (cm)');
ylabel('Piezometer head (cm)');
title('Hydraulic gradient line for observation 2')
grid on
subplot(2,2,3)
x=[0,6.028,6.868,7.318,8.108,14.154]
y=[22.501,22.0107,21.0206,20.0305,18.0401,20.00102]
plot(x,y)
xlabel('Distance of section from section a (cm)');
ylabel('Total head (cm)');
title('Total energy gradient line for observation 2')
grid on
subplot(2,2,4)
x=[0,6.028,6.868,7.318,8.108,14.154]
y=[25.001,24.0107,23.0206,22.0305,20.5401,22.30102]
plot(x,y)
xlabel('Distance of section from section a (cm)');
ylabel('Total head (cm)');
title('Total energy gradient line for observation 1')
grid on
hold off
