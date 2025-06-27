% RK4 method solver for vehicle dynamics - Grid independence study

clear; clc;

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
T = 5;             % Total simulation time (s)
dt = [0.1, 0.05, 0.02, 0.01, 0.005]; % Step Values

% Define A and B matrices (constant for constant u)
A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
     (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];

B = [Caf/m; a*Caf/Iz];

% Define ODE function: dx/dt = A*x + B*delta
f = @(t, x) A * x + B * delta;

% Initial condition: lateral velocity and yaw rate both zero
x0 = [0; 0];

% Initialize a struct array to store results for each iteration
results = struct();

% Loop over dt values
for i = 1:length(dt)
    h = dt(i); % Current step size
    t = 0:h:T; % Time vector for current step size
    
    % Solve IVP for current step size
    [t, x] = solveIVP(f, [0, T], x0, h, @rk4);
    
    % Store results in the struct
    results(i).dt = h; % Store the step size
    results(i).t = t;  % Store the time vector
    results(i).x = x;  % Store the solution matrix
    
    % Plot results for current step size
    figure;
    plot(t, x(1,:), 'b', 'LineWidth', 2);
    title(['Time vs Yaw for dt = ', num2str(h)]);
    ylabel('Lateral Velocity');
    xlabel('Time');
    grid on;
end

t = 0:1:T;
A_sol = [-13.0964*exp(-1.9745*t) + 24.4684*exp(-0.9839*t) - 11.3720; 
    -0.2496*exp(-1.9745*t)-0.6962*exp(-0.9839*t)+0.9457];

[L2LatSpeedNorm, L2YawRateNorm] = grid_check_L2_norm(f, [0,T], x0, dt, A_sol(1,:), A_sol(2,:))

% ==== FUNCTION DEFINITIONS ====

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

function [L2LatSpeedNorm, L2YawRateNorm] = grid_check_L2_norm(ODE, tspan, x0, stepSizes, actLatSpeed, actYawRates)
    L2LatSpeedNorm = zeros(length(stepSizes));
    L2YawRateNorm = zeros(length(stepSizes));
    times = tspan(1):1:tspan(2);

    for i = 1:length(stepSizes)
        [t, x_dot] = solveIVP(ODE, tspan, x0, stepSizes(i), @rk4);
        latSpeedValues = zeros(0, length(times)); %vector to store length(times) values
        yawRateValues = zeros(0, length(times));

        for j = 1:length(times)
            temp_indx = find(t == times(j));
            latSpeedValues(j) = x_dot(1,temp_indx)
            yawRateValues(j) = x_dot(2, temp_indx)
        end 
        L2_norm_lat = sqrt(sum(actLatSpeed - latSpeedValues).^2);
        L2_norm_yaw = sqrt(sum(actYawRates - yawRateValues).^2);
        
        L2LatSpeedNorm(i) = L2_norm_lat;
        L2YawRateNorm(i) = L2_norm_yaw;
    end
    length(L2YawRateNorm)
    figure (1)
    plot(log10(stepSizes), log10(L2LatSpeedNorm), 'o-');
    title('Lateral speed grid refinement')
    xlabel('log(step size)')
    ylabel('log(error)')

    figure (2)
    plot(log10(stepSizes), log10(L2YawRateNorm), 'o-');
    title('Yaw rate grid refinement');
    xlabel('log(step size)')
    ylabel('log(error)')
end
