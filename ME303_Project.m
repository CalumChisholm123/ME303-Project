

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
%RK4 Plots for Part A1
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


%=========Euler Method (forward)================
function [t,x] = euler_1(ode, tspan, x0, dt)

    t = tspan(1) : dt : tspan(2);
    n = length(t);
    x = zeros(length(x0), n);
    x(:,1) = x0;

    for i = 1:length(t) -1;
        x(:, i+1) = x(:, i) + dt * ode(t(i), x(:, i));
    end
end
%========Euler Plots for Part A1=================
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
ylabel('Yaw Rate (rad/s)')
title("Time Vs. Yaw Rate - Euler");

%B1 RK4 Plot
u_values = [20,50,75,100,200,300] / 3.6;

figure(5);
hold on;
for i = 1:length(u_values)
    u = u_values(i);

    A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
     (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];
    
    B = [Caf/m; a*Caf/Iz];

    f = @(t, x) A * x + B * delta;
    [t, x] = solveIVP(f, [0, T], x0, dt, @rk4);

    plot(t,x(1,:),'Displayname', ['u = ',num2str(u), 'm/s']);

end  
legend;
xlabel('Times (s)');
ylabel ('Lateral Acceleration (m/s2)');
title("Lateral Acceleration Vs. Various Speeds - RK4");
figure(6);
hold on;
for i = 1:length(u_values)
    u = u_values(i);

    A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
     (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];
    
    B = [Caf/m; a*Caf/Iz];

    f = @(t, x) A * x + B * delta;
    [t, x] = solveIVP(f, [0, T], x0, dt, @rk4);
   
    plot(t,x(2,:),'Displayname', ['u = ',num2str(u), 'm/s']);

end  
legend;
xlabel('Times (s)');
ylabel ('Yaw (rad/s)');
title("Yaw Rate Vs. Various Speeds - RK4");

%% ============ Part A2 ============
%Confirm the order of solvers
stepSizes = [0.1, 0.05, 0.01, 0.001];
[euler_errors, RK4_errors] = grid_check_L2_norm(f, tspan, x0, stepSizes, 0.0000001); %10^-7 as ground truth 
euler_slope = polyfit(log(stepSizes), log(euler_errors), 1);
RK4_slope = polyfit(log(stepSizes), log(RK4_errors), 1);

%Plot points at t = 1
[lateral_vels, yaw_rates] = grid_residuals(f, tspan, x0, stepSizes, 1);

function [euler_errors, RK4_errors] = grid_check_L2_norm(ODE, tspan, x0, stepSizes, truth_size)
    euler_errors = zeros(length(stepSizes),1); %Store L2 values for different stepSizes
    RK4_errors = zeros(length(stepSizes),1);

    %Solve the ground truth for t = 1 s and RK4 method
    time_select = 1;
    [t_true, x_true] = solveIVP(ODE, tspan, x0, truth_size, @rk4);
    t_index = find(t_true==time_select);
    tru_lat_speed = x_true(1, t_index);
    tru_yaw_rate = x_true(2, t_index);

    for i = 1:length(stepSizes) %iterate through each step size (h) value
        [t, x_dot] = euler_1(ODE, tspan, x0, stepSizes(i));
        temp_index = find(t == time_select);
        lat_speed = x_dot(1, temp_index);
        yaw_rate = x_dot(2, temp_index);
        euler_errors(i) = sqrt((tru_lat_speed - lat_speed)^2 + (tru_yaw_rate - yaw_rate)^2);

        %Do the same for RK4 error
        [t, x_dot] = solveIVP(ODE, tspan, x0, stepSizes(i), @rk4);
        temp_index = find(t == time_select);
        lat_speed = x_dot(1, temp_index);
        yaw_rate = x_dot(2, temp_index);
        RK4_errors(i) = sqrt((tru_lat_speed - lat_speed)^2 + (tru_yaw_rate - yaw_rate)^2);
    end
    
    figure (1)
    plot(log10(stepSizes), log10(euler_errors),'DisplayName',"Euler's");
    hold on
    plot(log10(stepSizes), log10(RK4_errors), 'DisplayName', 'RK4');
    title("Euler's and RK4 Error")
    xlabel('log(grid size)')
    ylabel('log(L2-norm error)')
    legend
    grid on
end

%Plot grid independence
function [lateral_vels, yaw_rates] = grid_residuals(ODE, tspan, x0, stepSizes, time_select)
    lateral_vels = zeros(2, length(stepSizes)); %1 for Euler, 2 for RK4
    yaw_rates = zeros(2, length(stepSizes));
    
    for i = 1:length(stepSizes)
        [t, x] = euler_1(ODE, tspan, x0, stepSizes(i));
        temp_indx = find(t == time_select);
        lateral_vels(1,i) = x(1,temp_indx);
        yaw_rates(1, i) = x(2, temp_indx);

        [t, x] = solveIVP(ODE, tspan, x0, stepSizes(i), @rk4);
        temp_index = find(t == time_select);
        lateral_vels(2, i) = x(1, temp_index);
        yaw_rates(2, i) = x(2, temp_index);
    end
end
%%
%============ Part B2 ============
clear all;
clc

% Parameters
m = 1400;          % kg
a = 0.988;          % distance of CM to front axle 
b = 1.482;          % distance of CM to rear axle 
Caf = 25000;       % N/rad
Car = 21000;       % N/rad
Iz = 2420;         % kg·m^2
delta = 0.1;       % Step steering input (rad)


u_values = (1:0.0001:300);

for i = 1:length(u_values)
    u = u_values(i);
    
    A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
     (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];
    
    lambda = eig(A);


    if any(real(lambda) >= 0)
        u = u*3.6;
        disp(u)

        disp('UNSTABLE: At least one eigenvalue has a non-negative real part.');

         break;
    end

 end



%============ Part B 3. ===========
u = 100/3.6; %chaning U to 100km/hr
tspan = [0 5]; %changing span to 5
L = a+b;

dt = 0.001; 
T_total = 50;%time that car is moving

figure(7);
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
hold on;
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

figure(8);
plot(u_values*3.6, ideal_list, 'k--', 'LineWidth', 2); hold on;
plot(u_values*3.6, actual_list, 'ro-', 'LineWidth', 2);
legend('Ideal Yaw Rate', 'Actual Yaw Rate');
xlabel('Speed (km/h)');
ylabel('Yaw Rate (rad)');
title('Ideal vs Actual Yaw Rate at Different Speeds');
grid on;

%============ Part C1 ====================
A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
     (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];

lambda = eig(A);
disp(lambda)

%% ==================== Section D1 ====================
Car_span = [10000, 50000];

[C_rears, stability_list, lambdas_real1, lambdas_real2, lambdas_imag1, lambdas_imag2] = check_stability(Caf, Car_span, 1, m, u, a, b, Iz);
Cars_list = [18000 21000 30000];
tspan = [0, 5];
[times, lat_vel_results, yaw_rate_results] = solve_diff_Cars(Cars_list, tspan, 0.001, x0, m, u, a, b, Iz, Caf, delta);

%Find lambda values for all rear cornering stiffness values
function [C_rears, stability_list, lambdas_real1, lambdas_real2, lambdas_imag1, lambdas_imag2] = check_stability(Caf, Car_span, stepSize, m, u, a, b, Iz)
    C_rears = Car_span(1) : stepSize : Car_span(2);
    stability_list = zeros(1, length(C_rears)); %0 of unstable, 1 if stable
    lambdas_real1 = zeros(1, length(C_rears));
    lambdas_real2 = zeros(1, length(C_rears));
    lambdas_imag1 = zeros(1, length(C_rears));
    lambdas_imag2 = zeros(1, length(C_rears));

    for i = 1:1:length(C_rears)
        Car = C_rears(i);
        A = [-(Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
              (-a*Caf + b*Car)/(Iz*u), - (a^2 * Caf + b^2 * Car)/(Iz*u)];

        lambdas = eig(A);
        lambdas_real1(i) = real(lambdas(1));
        lambdas_real2(i) = real(lambdas(2));
        lambdas_imag1(i) = imag(lambdas(1));
        lambdas_imag2(i) = imag(lambdas(2));
    end
end

%Solve lateral velocity and yaw rates for different rear cornering
%stiffnesses and plot them
function [times, lat_vel_results, yaw_rate_results] = solve_diff_Cars(Cars_list, tspan, stepSize, x0, m, u, a, b, Iz, Caf, delta)
    times = tspan(1):stepSize:tspan(2);
    lat_vel_results = zeros(length(Cars_list), length(times));
    yaw_rate_results = zeros(length(Cars_list), length(times));

    for i = 1:length(Cars_list)
        Car = Cars_list(i);
        A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
            (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];
        B = [Caf/m; a*Caf/Iz];

        f = @(t, x) A * x + B * delta;

        [t, x] = solveIVP(f, tspan, x0, stepSize, @rk4);
        lat_vel_results(i, :) = x(1, :);
        yaw_rate_results(i, :) = x(2, :);
    end

    %Plot everything in a separate loop for clarity
    figure(1)
    for i = 1:length(Cars_list)
        plot(times, lat_vel_results(i,:), "DisplayName", sprintf("C_{\\alphar} = %.0f N/rad", Cars_list(i)));
        hold on
    end
    title('$\dot{y}$ (m/s) for constant $C_{\alpha f}$ and varying $C_{\alpha r}$', 'Interpreter', 'latex');
    ylabel('$\dot{y}$ (m/s)', 'Interpreter', 'latex')
    xlabel("Time (s)")
    legend
    grid on

    figure(2)
    for i = 1:length(Cars_list)
        plot(times, yaw_rate_results(i,:), "DisplayName", sprintf("C_{\\alphar} = %.0f N/rad", Cars_list(i)));
        hold on
    end
    title('$\dot{\psi}$ (rad/s) for constant $C_{\alpha f}$ and varying $C_{\alpha r}$', 'Interpreter', 'latex');
    ylabel('$\dot{\psi}$ (m/s)', 'Interpreter', 'latex')
    xlabel("Time (s)")
    legend
    grid on
  end

%% ==================== Section D2 ====================

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

%============ Part F =================
m = 1580;          % kg
a = 1.1836;          % m
b = 1.5064;          % m
Caf = 69800;       % N/rad
Car = 69900;       % N/rad
Iz = 2817.1;         % kg·m^2
u = 180 / 3.6;      % Convert from km/h to m/s
delta = 0.1;       % Step steering input (rad)

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
title('Handling Behaviour of RAV4');
grid on;