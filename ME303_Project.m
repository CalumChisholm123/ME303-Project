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

%% ==================== Section B.1 ====================

% Calums Speed Changer 
u_values = [20,50,75,100,200,300] / 3.6;

%RK4 different Longitudinal Speeds
for i = 1:length(u_values)
    u = u_values(i);

    A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
     (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];

    B = [Caf/m; a*Caf/Iz];

    f = @(t, x) A * x + B * delta;
    [t, x] = solveIVP(f, [0, T], x0, dt, @rk4);

    figure(2);
    hold on;
    plot(t,x(1,:),'DisplayName', ['u = ',num2str(u), 'm/s']);
    legend('Location','best');
    xlabel('Times (s)');
    ylabel('Lateral Acceleration (m/s)')
    title("Lateral Acceleration vs Time for Different u values - RK4");

    figure(3);
    hold on;
    plot(t,x(2,:),'DisplayName', ['u = ',num2str(u), 'm/s']);
    legend('Location','best');
    xlabel('Times (s)');
    ylabel('Yaw Rate')
    title("Yaw Rate vs Time for Different u values - RK4");

end  

%Eulers Longitudinal Speeds
for p = 1:length(u_values)
    u = u_values(p);

    A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
     (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];

    B = [Caf/m; a*Caf/Iz];

    f = @(t, x) A * x + B * delta;
    [t, x] = euler_1(f, tspan, x0, dt);

    figure(4);
    hold on;
    plot(t,x(1,:),'DisplayName', ['u = ',num2str(u), 'm/s']);
    legend;
    xlabel('Times (s)');
    ylabel('Lateral Acceleration (m/s2)')
    title("Lateral Acceleration vs Time for Different u values - Euler");

    figure(5);
    hold on;
    plot(t,x(2,:),'DisplayName', ['u = ',num2str(u), 'm/s']);
    legend;
    xlabel('Times (s)');
    ylabel('Yaw Rate')
    title("Yaw Rate vs Time for Different u values - Euler");

 end  

%% ==================== Section B.2 ====================

T = 10000; % Total simulation time (s)
dt = 0.1; % Step Values 
N = T/dt; % Number of Steps 
t = linspace(0,T,N+1); % Time Vector
u = 240 / 3.6 % Convert from km/h to m/s
threshold = 999999; %arbitrary large slope   
  
for u = 200/3.6:250/3.6

    A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
    (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];

    B = [Caf/m; a*Caf/Iz];

    f = @(t, x) A * x + B * delta;
    [t, x] = solveIVP(f, [0, T], x0, dt, @rk4);


    dx_dt = diff(x, 1, 2) / dt;  % Approximate time derivative of x
    max_slope = max(abs(dx_dt), [], 'all');  % Max absolute slope across both states

    if max_slope > threshold
        warning(['Diverging slope detected at u = ', num2str(u)]);
        break;
    end

    figure(6);
    hold on;
    plot(t,x(1,:),'DisplayName', ['u = ',num2str(u), 'm/s']);
    legend('Location','best');
    xlabel('Times (s)');
    ylabel('Lateral Acceleration (m/s)')
    title("Lateral Acceleration vs Time for Different u values - RK4");

    figure(7);
    hold on;
    plot(t,x(2,:),'DisplayName', ['u = ',num2str(u), 'm/s']);
    legend('Location','best');
    xlabel('Times (s)');
    ylabel('Yaw Rate')
    title("Yaw Rate vs Time for Different u values - RK4");

end
   
%% ==================== Section B.3 ====================

T = 10000; % Total simulation time (s)
dt = 0.01; % Step Values 
N = T/dt; % Number of Steps 
t = linspace(0,T,N+1); % Time Vector
u = 229 / 3.6; % Convert from km/h to m/s

A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
    (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];
B = [Caf/m; a*Caf/Iz];
f = @(t, x) A * x + B * delta;

[t_track, x_track] = solveIVP(f, [0, T], x0, dt, @rk4);

% Extract lateral velocity and yaw rate
v = x_track(1,:);  % lateral velocity (m/s)
r = x_track(2,:);  % yaw rate (rad/s)

% Integrate yaw angle over time to get heading direction
psi = cumtrapz(t, r);  % heading angle (rad)

% Compute global velocities in X and Y directions
x_dot = u * cos(psi) - v .* sin(psi);
y_dot = u * sin(psi) + v .* cos(psi);

% Integrate to get global positions
X = cumtrapz(t, x_dot);
Y = cumtrapz(t, y_dot);

% Plot vehicle trajectory (X vs Y)
figure(8);
plot(X, Y);
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('Vehicle Track');
grid on;
axis equal;

%% ==================== Section C2 ====================
clear all; 
close all;
clc;

m = 1400;          % kg
a = 1.14;          % m
b = 1.33;          % m
Caf = 5000;       % N/rad
Car = 5000;       % N/rad
Iz = 2420;         % kg·m^2
u = 100 / 3.6;      % Convert from km/h to m/s
delta = 0.1;       % Step steering input (rad)

% Time setup
tspan = [0,30];
T = 60;             % Total simulation time (s)
dt = 0.01; % Step Values 
N = T/dt; % Number of Steps 
t = linspace(0,T,N+1); % Time Vector 

A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
    (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];
B = [Caf/m; a*Caf/Iz];

f = @(t, x) A * x + B * delta;
x0 = [0; 0];


[t_track, x_track] = solveIVP(f, [0, T], x0, dt, @rk4);

% Extract lateral velocity and yaw rate
v = x_track(1,:);  % lateral velocity (m/s)
r = x_track(2,:);  % yaw rate (rad/s)

% Integrate yaw angle over time to get heading direction
psi = cumtrapz(t, r);  % heading angle (rad)

% Compute global velocities in X and Y directions
x_dot = u * cos(psi) - v .* sin(psi);
y_dot = u * sin(psi) + v .* cos(psi);

% Integrate to get global positions
X = cumtrapz(t, x_dot);
Y = cumtrapz(t, y_dot);

% Plot vehicle trajectory (X vs Y)
figure(9);
plot(X, Y);
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('Vehicle Track - Normal Conditions - Normal Tires');
grid on;
axis equal;