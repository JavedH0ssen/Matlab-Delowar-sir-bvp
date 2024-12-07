
% Plot the solution
plot(sol.x, sol.y(1,:));
hold on;
plot(sol.x, sol.y(2,:));
xlabel('x');
ylabel('y');
legend('y1', 'y2');

% Solve the BVP
sol = bvp4c(@dydx, @bc, @solinit, [0, pi]);
function dydx = dydx(x,y)
    dydx = [y(2); -y(1) - sin(x)];
end

function bc = bc(ya,yb)
    bc = [ya(1); ya(2) - 1; yb(1); yb(2)];
end

function solinit = solinit(x)
    solinit = [0; 0; 1; 0];
end


