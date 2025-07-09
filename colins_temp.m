clear all; 
close all;
clc;

% Parameters
m = 1400;          % kg
a = 1.14;          % m
b = 1.33;          % m
Caf = 25000;       % N/rad
Car = 21000;       % N/rad
Iz = 2420;         % kg·m^2
u = 75 / 3.6;      % Convert from km/h to m/s
delta = 0.1;       % Step steering input (rad)

% Time setup
tspan = [0,30];
T = 5;             % Total simulation time (s)
dt = 0.01; % Step Values 
N = T/dt; % Number of Steps 
t = linspace(0,T,N+1); % Time Vector 



% Define A and B matrices (constant for constant u)
A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
     (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];

B = [Caf/m; a*Caf/Iz];

% Define ODE function: dx/dt = A*x + B*delta
f = @(t, x) A * x + B * delta;

% Initial condition: lateral velocity and yaw rate both zero
x0 = [0; 0];

% ==== FUNCTION DEFINITIONS ====

[t, x] = solveIVP(f, [0, T], x0, dt, @rk4);



%============ Part B 3. ===========
u = 75/3.6; %changing from km/h to m/s
L = a+b;

x0 = [0; 0];
% Diff steering inputs
delta_values = [0.1];

dt = 0.01; 
T_total = 10;%time that car is moving

figure;
hold on;

for i=1:length(delta_values)

    delta = delta_values(i);
    A = [-(Caf + Car)/(m*u), (a*Caf + b*Car)/(m*u) - u;
         (-a*Caf + b*Car)/(Iz*u), -(a^2*Caf + b^2*Car)/(Iz*u)];

    B = [Caf/m;
         a*Caf/Iz];

    f = @(t, x) A * x + B * delta;
    [t_track, x_track] = solveIVP(f, [0, T_total], x0, dt, @rk4);

    v = x_track(1,:);
    r = x_track(2,:);

    phi = zeros(1, length(t_track));%Riemann Sum
    for i = 2:length(t_track)
        phi(i) = phi(i-1) + r(i-1) * (t_track(i) - t_track(i-1));
    end

    x_dot = u * cos(phi) - (v + a.*r) .* sin(phi);
    y_dot = u * sin(phi) + (v + a.*r) .* cos(phi);

    X = zeros(1, length(t_track));%Integration
    Y = zeros(1, length(t_track));
    for i = 2:length(t_track)
        X(i) = X(i-1) + x_dot(i-1) * (t_track(i) - t_track(i-1));
        Y(i) = Y(i-1) + y_dot(i-1) * (t_track(i) - t_track(i-1));
    end

    plot(X, Y, 'LineWidth', 2, 'DisplayName', sprintf('%.1f rad', delta));
    hold on
end

xlabel('X (m)');
ylabel('Y (m)');
title('Handling Behaviour of Car with Step Steering Experiment');
legend;
grid on;
axis equal;

% Generic IVP solver
function [t, y] = solveIVP(f, tspan, y0, h, solver)
    t = tspan(1):h:tspan(2);
    y = zeros(length(y0), length(t));
    y(:,1) = y0;
    for n = 1:length(t) - 1
        y(:,n+1) = solver(f, t(n), y(:,n), h);
    end
end

% Runge-Kutta 4th order method
function ynew = rk4(f, t, y, h)
    k1 = f(t, y);
    k2 = f(t + 0.5 * h, y + 0.5 * h * k1);
    k3 = f(t + 0.5 * h, y + 0.5 * h * k2);
    k4 = f(t + h, y + h * k3);
    ynew = y + (h / 6) * (k1 + 2*k2 + 2*k3 + k4);
end

%Euler Method (forward)
function [t,x] = euler_1(ode, tspan, x0, dt)

    t = tspan(1) : dt : tspan(2);
    n = length(t);
    x = zeros(length(x0), n);
    x(:,1) = x0;

    for i = 1:length(t) -1;
        x(:, i+1) = x(:, i) + dt * ode(t(i), x(:, i));
    end
end