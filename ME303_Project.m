% RK4 method solver for vehicle dynamics - Grid independence study

clear; clc;

% Parameters
m = 1400;          % kg
a = 1.14;          % m
b = 1.33;          % m
Caf = 25000;       % N/rad
Car = 21000;       % N/rad
Iz = 2420;         % kg·m^2
u = 330 / 3.6;      % Convert from km/h to m/s
delta = 0.1;       % Step steering input (rad)

% Time setup
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

plot(t,x(1,:),"b", LineWidth=2);
title("Time vs Yaw")
ylabel("Lateral Velocity");
xlabel("Time");
grid on;


%lateralspeed = -13.0964*exp(-1.9745*t) + 24.4684*exp(-0.9839*t) - 11.3720;
% yaw = -0.2496*exp(-1.9745*t) - 0.6962*exp(-0.9839*t) + 0.9457;
% plot(t,lateralspeed)
% title("Time vs Lateral Speed Sol")
% plot(t,yaw)
% title("Time vs Yaw Solution")