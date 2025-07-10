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

%% ==================== FUNCTION DEFINITIONS ====================

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

figure(1);
hold on;
plot(t,x(1,:),"b", LineWidth=2);
title("Time vs Yaw")
ylabel("Lateral Velocity");
xlabel("Time");
grid on;


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


%% ==================== Section D2 ====================

clear all;
close all;
clc;

% Vehicle Parameters
a = 1.14;          % m
b = 1.33;          % m
Caf = 5000;       % N/rad
Car = 5000;       % N/rad
Iz = 2420;         % kg·m^2
u = 100 / 3.6;     % Convert from km/h to m/s
delta = 0.1;       % Step steering input (rad)
m_list = 700:100:1400; % List of masses to simulate (kg)

% Time setup
tspan = [0, 60];   % Total simulation time (s)
dt = 0.001;         % Time step
t = tspan(1):dt:tspan(2); % Time vector

% Create figure and hold for overlaid plots
figure(1);
hold on;

% Loop through each mass in the list
for i = 1:length(m_list)
    m = m_list(i); % Current mass

    A = [-(Caf + Car)/(m*u), (a*Caf + b*Car)/(m*u) - u;
    (-a*Caf + b*Car)/(Iz*u), -(a^2*Caf + b^2*Car)/(Iz*u)];

    B = [Caf/m;
        a*Caf/Iz];
    
    f = @(t, x) A * x + B * delta;

    x0 = [0; 0];

    [t_track, x_track] = solveIVP(f, [0, 120], x0, dt, @rk4);
    
    v = x_track(1,:);
    r = x_track(2,:);
    
    psy = zeros(1, length(t_track));%Riemann Sum
    for i = 2:length(t_track)
        psy(i) = psy(i-1) + r(i-1) * (t_track(i) - t_track(i-1));
    end
    
    x_dot = u * cos(psy) - (v + a.*r) .* sin(psy);
    y_dot = u * sin(psy) + (v + a.*r) .* cos(psy);
    
    X = zeros(1, length(t_track));%Integration
    Y = zeros(1, length(t_track));
    for i = 2:length(t_track)
        X(i) = X(i-1) + x_dot(i-1) * (t_track(i) - t_track(i-1));
        Y(i) = Y(i-1) + y_dot(i-1) * (t_track(i) - t_track(i-1));
    end

    plot(X, Y, 'DisplayName', sprintf('Mass = %d kg', m));
end

% Finalize the plot
hold off;
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('Vehicle Track for Different Masses');
grid on;
axis equal;
legend('Location', 'northwest');
