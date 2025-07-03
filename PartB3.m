

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
    figure(1);
    hold on;
    plot(t,x(1,:),'DisplayName', ['u = ',num2str(u), 'm/s']);
    legend;
    xlabel('Times (s)');
    ylabel('Lateral Acceleration (m/s)')
    title("Different u values - RK4");
    figure(2);
    hold on;
    plot(t,x(2,:),'DisplayName', ['u = ',num2str(u), 'm/s']);
    legend;
    xlabel('Times (s)');
    ylabel('Yaw Rate')
    title("Different u values - RK4");

end  

%Eulers Longitudinal Speeds
for p = 1:length(u_values)
    u = u_values(p);

    A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
     (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];
    
    B = [Caf/m; a*Caf/Iz];

    f = @(t, x) A * x + B * delta;

    [t, x] = euler_1(f, tspan, x0, dt);

    figure(3);
    hold on;
    plot(t,x(1,:),'DisplayName', ['u = ',num2str(u), 'm/s']);
    legend;
    xlabel('Times (s)');
    ylabel('Lateral Acceleration (m/s2)')
    title("Different u values - Euler");
    figure(4);
    hold on;
    plot(t,x(2,:),'DisplayName', ['u = ',num2str(u), 'm/s']);
    legend;
    xlabel('Times (s)');
    ylabel('Yaw Rate')
    title("Different u values - Euler");

end  
%============ Part B 3. ===========
u = 100; %chaning U to 100km/hr
tspan = [0 5]; %changing span to 5


[t_track , x_track] = solveIVP(f, tspan, x0, dt, @rk4);
v = x_track(1,:);
r = x_track(2,:);

direction = cumtrapz(t_track, r);

x_dot = u * cos(direction) - (v + (r*a)) .* sin(direction);
y_dot = u * sin(direction) + (v + (r*a)) .* cos(direction);

X = cumtrapz(t_track, x_dot);
Y = cumtrapz(t_track, y_dot);


figure(5);
plot(X,Y, 'LineWidth',2);
xlabel ('X Position (m)');
ylabel ('Y position (m)');
title ('Vehicle Track at 100km/hr');
grid on;
axis equal;

