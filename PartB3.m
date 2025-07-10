

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

figure(1)
plot(t,x(2,:),"b", LineWidth=2);
hold on;
title("Time vs Yaw - RK4")
ylabel("Yaw Rate (degrees)");
xlabel("Time");
grid on;

figure(2)
plot(t,x(1,:),"b", LineWidth=2);
hold on;
title("Time vs Lateral Velocity - RK4")
ylabel("Lateral Velocity (m/s)");
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

figure(3);
hold on;
plot(t,x(1,:), "b", LineWidth=2);
xlabel('Times (s)');
ylabel('Lateral Acceleration (m/s2)')
title("Time Vs. Lateral Acceleration - Euler");
figure(4);
hold on;
plot(t,x(2,:), "b", LineWidth=2);
xlabel('Times (s)');
ylabel('Yaw Rate')
title("Time Vs. Yaw Rate - Euler");
% Calums Speed Changer 
% u_values = [20,50,75,100,200,300] / 3.6;

% %RK4 different Longitudinal Speeds
% for i = 1:length(u_values)
%     u = u_values(i);
% 
%     A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
%      (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];
% 
%     B = [Caf/m; a*Caf/Iz];
% 
%     f = @(t, x) A * x + B * delta;
%     [t, x] = solveIVP(f, [0, T], x0, dt, @rk4);
%     figure(1);
%     hold on;
%     plot(t,x(1,:),'DisplayName', ['u = ',num2str(u), 'm/s']);
%     legend;
%     xlabel('Times (s)');
%     ylabel('Lateral Acceleration (m/s)')
%     title("Different u values - RK4");
%     figure(2);
%     hold on;
%     plot(t,x(2,:),'DisplayName', ['u = ',num2str(u), 'm/s']);
%     legend;
%     xlabel('Times (s)');
%     ylabel('Yaw Rate')
%     title("Different u values - RK4");

% end  





% end  

%============ Part B 3. ===========
u = 100/3.6; %chaning U to 100km/hr
tspan = [0 5]; %changing span to 5
L = a+b;
u_values = [20, 40, 60, 80, 100, 200] / 3.6;

yaw_rate = 0.1:0.1:0.5;

x0 = [0; 0];


dt = 0.001; 
T_total = 50;%time that car is moving

figure;
hold on;



A = [-(Caf + Car)/(m*u), (a*Caf + b*Car)/(m*u) - u;
    (-a*Caf + b*Car)/(Iz*u), -(a^2*Caf + b^2*Car)/(Iz*u)];

B = [Caf/m;
    a*Caf/Iz];

f = @(t, x) A * x + B * delta;
[t_track, x_track] = solveIVP(f, [0, T_total], x0, dt, @rk4);

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

plot(X, Y, 'LineWidth', 2, 'DisplayName', sprintf('%.1f rad', delta));
hold on


xlabel('X (m)');
ylabel('Y (m)');
title('Handling Behaviour of Car with Step Steering Experiment');
grid on;



ideal_list = [];
actual_list = [];

for i = 1:length(u_values)
    u = u_values(i);

    A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
         (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];
    B = [Caf/m; a*Caf/Iz];

    f = @(t, x) A * x + B * delta;
    x0 = [0; 0];

    T = 10;
    dt = 0.01;
    [t, x] = solveIVP(f, [0 T], x0, dt, @rk4);

    ideal = (u * delta) / L;
    actual = mean(x(2,end-100:end)); % average of last ~1s to get steady-state

    ideal_list(end+1) = ideal;
    actual_list(end+1) = actual;
end

figure;
plot(u_values*3.6, ideal_list, 'k--', 'LineWidth', 2); hold on;
plot(u_values*3.6, actual_list, 'ro-', 'LineWidth', 2);
legend('Ideal Yaw Rate', 'Actual Yaw Rate');
xlabel('Speed (km/h)');
ylabel('Yaw Rate (rad)');
title('Ideal vs Actual Yaw Rate at Different Speeds');
grid on;