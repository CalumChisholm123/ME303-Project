%Parameters for Vehicle
m = 1400;          % kg
a = 1.14;          % m
b = 1.33;          % m
Caf = 25000;       % N/rad
Car = 21000;       % N/rad
Iz = 2420;         % kg·m^2
u = 75 / 3.6;      % Convert from km/h to m/s
delta = 0.1;       % Step steering input (rad)

%Matrix
A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
     (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];

B = [Caf/m;a*Caf/Iz];

bicycle_ODE = @(t,x) A * x + B * delta;

%Define IVP
tspan = [0,5];
y0 = 1;
x0 = [0,0];
h = 0.01;
f = @(t,y) t * y;

% Generic IVP solver
function [t, y] = solveIVP(f, tspan, y0, h, solver)
    t = tspan(1):h:tspan(2);
    y = zeros(length(y0), length(t));
    y(:,1) = y0;
    for n = 1:length(t) - 1
        y(:,n+1) = solver(f, t(n), y(:,n), h);
    end
end

%Euler Method
function [t,x] = euler_1(ode, tspan, x0, h)

    t = tspan(1) : h : tspan(2);
    n = length(t);
    x = zeros(length(x0), n);
    x(:,1) = x0;

    for i = 1:length(t) -1;
        x(:, i+1) = x(:, i) + h * ode(t(i), x(:, i));
    end
end

stepSize = [0.1, 0.05, 0.01, 0.001];

figure;
hold on;
for h = stepSize
    [t, x] = euler_1(bicycle_ODE, tspan, x0, h);
    plot(t, x(2,:), 'DisplayName', ['h = ', num2str(h)]);
end

title('Euler Method - Grid Independence Check');
xlabel('Time (s)');
ylabel('Yaw rate (rad/s)');
legend;
grid on;

h_selected = 0.01; % Chosen based on grid study
[t, x] = euler_method(bicycle_ode, tspan, x0, h_selected);

figure;
plot(t, x(2,:), 'b', 'LineWidth', 1.5);
title('Final Euler Solution (h = 0.01)');
xlabel('Time (s)');
ylabel('Yaw rate (rad/s)');
grid on;