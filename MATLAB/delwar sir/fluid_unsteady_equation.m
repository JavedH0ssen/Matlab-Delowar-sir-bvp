clear
tic;
iteration=1;
n_points=51;
dom_size=1;
h=dom_size/(n_points-1);
dt=0.00001;
alpha=dt/(h*h);

y(n_points,n_points)=0;
y(1,:)=1;
y(:,1)=1;
y(:,n_points)=1;
y(n_points,:)=5;
ynew=y;
error_mag=1;N=0;
error_required=1e-3;
y_transient_C=cell(1,1);


while error_mag>error_required
    yold = ynew;  % Store the old values of y
    % Vectorized computation of ynew
    ynew(2:n_points-1, 2:n_points-1) = y(2:n_points-1, 2:n_points-1) + ...
        alpha * (y(3:n_points, 2:n_points-1) + y(1:n_points-2, 2:n_points-1) + ...
        y(2:n_points-1, 3:n_points) + y(2:n_points-1, 1:n_points-2) - ...
        4*y(2:n_points-1, 2:n_points-1));
    N = N + (n_points-2)^2;  % Update N
    % Vectorized computation of error_mag
    error_mag = sum(sum(abs(ynew(2:n_points-1, 2:n_points-1) - yold(2:n_points-1, 2:n_points-1))));
    y_transient_C{iteration}=ynew;  
    iteration=iteration+1;
    y=ynew;
end

N
iteration
tyme=toc;
fprintf('\nElapsed time is %.6f seconds.\n', tyme);
%% plotting
% timeSelected=100;
x_dom=((1:n_points)-1)*h;
y_dom=(1-((1:n_points)-1))*h;
[X,Y]=meshgrid(x_dom,y_dom);
% ytimestep=y_transient(timeSelected,:,:);
% ytimestep=reshape(ytimestep,[n_points,n_points]);
% contourf(X,Y,ytimestep)
% colorbar
%% Animation after Every N steps
Nn=20;
figure;
pause(02)
for i=1:Nn:iteration
    ytimestep=cell2mat(y_transient_C(i));
    ytimestep=reshape(ytimestep,[n_points,n_points]);
    contourf(X,Y,ytimestep)
    colorbar
    title(['Time Step: ',num2str(i*dt)])
    pause(0.00001)
end
