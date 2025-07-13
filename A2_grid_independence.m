%Parameters for Vehicle
m = 1400;          % kg
a = 1.14;          % m
b = 1.33;          % m
Caf = 25000;       % N/rad
Car = 21000;       % N/rad
Iz = 2420;         % kg·m^2
u = 75 / 3.6;      % Convert from km/h to m/s
delta = 0.1;       % Step steering input (rad)

% Define A and B matrices (constant for constant u)
A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
     (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];

B = [Caf/m; a*Caf/Iz];

f = @(t, x) A * x + B * delta;
x0 = [0; 0];

%Define IVP
tspan = [0,5];
y0 = 1;
x0 = [0 0];
stepSizes = [0.1, 0.05, 0.01, 0.001];

%Confirm the order of solvers
[euler_errors, RK4_errors] = grid_check_L2_norm(f, tspan, x0, stepSizes, 0.0000001); %10^-7 as ground truth 
euler_slope = polyfit(log(stepSizes), log(euler_errors), 1);
RK4_slope = polyfit(log(stepSizes), log(RK4_errors), 1);

%Plot points at t = 1
[lateral_vels, yaw_rates] = grid_residuals(f, tspan, x0, stepSizes, 1);

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

    %Plot out the error
    %Lateral velocity 
    figure(1)
    plot(stepSizes, lateral_vels(1,:), "DisplayName", "Euler's")
    hold on
    plot(stepSizes, lateral_vels(2, :), "DisplayName", "RK4")
    xlabel("Step size")
    ylabel("Lateral velocity (m/s)")
    title("Lateral velocity vs step size at time = 1s")
    legend
    grid on

    figure(2)
    plot(stepSizes, yaw_rates(1, :), "DisplayName", "Euler's")
    hold on
    plot(stepSizes, yaw_rates(2, :), "DisplayName", "RK4")
    xlabel("Step size")
    ylabel("Yaw rate (rad/s)")
    title("Yaw rate vs step size at time = 1s")
    legend 
    grid on
end

