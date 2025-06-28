% RK4 method solver for vehicle dynamics - Verification with Exact Solution

clear; clc; close all;

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
step_sizes = [0.1, 0.05, 0.02, 0.01, 0.005]; % Step Values

% Define A and B matrices (constant for constant u)
A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
     (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];
B = [Caf/m; a*Caf/Iz];

% Define ODE function: dx/dt = A*x + B*delta
f = @(t, x) A * x + B * delta;

% Define the exact analytical solution as a function of time
f_exact = @(t) [-13.0964*exp(-1.9745*t) + 24.4684*exp(-0.9839*t) - 11.3720; 
                -0.2496*exp(-1.9745*t) - 0.6962*exp(-0.9839*t) + 0.9457];

% Initial condition: lateral velocity and yaw rate both zero
x0 = [0; 0];

% ==============================================================
% ==== SIMULATION AND ERROR CALCULATION ========================
% ==============================================================

% Initialize a struct array to store results and an array for errors
results = struct();
errors = zeros(size(step_sizes));

% Calculate the true solution at the final time T for error comparison
x_true_final = f_exact(T);

% Loop over dt values
fprintf('Running simulations and calculating errors...\n');
for i = 1:length(step_sizes)
    h = step_sizes(i); % Current step size
    
    % Solve IVP for current step size
    [t, x] = solveIVP(f, [0, T], x0, h, @rk4);
    
    % Store results in the struct
    results(i).dt = h; % Store the step size
    results(i).t = t;  % Store the time vector
    results(i).x = x;  % Store the solution matrix
    
    % Calculate the true error at T=5s using the Euclidean norm
    errors(i) = norm(x(:,end) - x_true_final);
end
fprintf('All simulations complete.\n\n');


% ==============================================================
% ==== CONVERGENCE ANALYSIS AND PLOTTING =======================
% ==============================================================

% 1. Calculate the observed order of accuracy 'p'
p = zeros(length(step_sizes)-1, 1);
for i = 1:length(errors)-1
    p(i) = log(errors(i) / errors(i+1)) / log(step_sizes(i) / step_sizes(i+1));
end

% 2. Display results in a table
fprintf('--- Verification Results ---\n');
fprintf('Step Size (h) |  True Error  | Observed Order (p)\n');
fprintf('---------------------------------------------------\n');
fprintf('%13.4f | %12.4e | \n', step_sizes(1), errors(1));
for i = 2:length(step_sizes)
    fprintf('%13.4f | %12.4e | %18.4f\n', step_sizes(i), errors(i), p(i-1));
end

% 3. Create a plot comparing numerical solutions to the exact solution
figure;
hold on;
% Plot the exact solution (Yaw Rate)
t_fine = 0:0.001:T;
x_exact_fine = f_exact(t_fine);
plot(t_fine, x_exact_fine(2,:), 'k-', 'LineWidth', 3, 'DisplayName', 'Exact Solution');

% Plot the numerical solutions
for i = 1:length(results)
    plot(results(i).t, results(i).x(2,:), '--', 'LineWidth', 1.5, ...
        'DisplayName', ['dt = ', num2str(results(i).dt)]);
end
hold off;
title('RK4 Grid Independence Study: Yaw Rate');
xlabel('Time (s)');
ylabel('Yaw Rate (rad/s)');
legend('show', 'Location', 'southeast');
grid on;

figure;
hold on;
% Plot the exact solution (Yaw Rate)
t_fine = 0:0.001:T;
x_exact_fine = f_exact(t_fine);
plot(t_fine, x_exact_fine(2,:), 'k-', 'LineWidth', 3, 'DisplayName', 'Exact Solution');

% Plot the numerical solutions
for i = 1:length(results)
    plot(results(i).t, results(i).x(2,:), '--', 'LineWidth', 1.5, ...
        'DisplayName', ['dt = ', num2str(results(i).dt)]);
end
hold off;
title('RK4 Grid Independence Study: Yaw Rate');
xlabel('Time (s)');
ylabel('Yaw Rate (rad/s)');
legend('show', 'Location', 'southeast');
grid on;



% Plot the exact solution (Lateral Speed)
t_fine = 0:0.001:T;
x_exact_fine = f_exact(t_fine);
plot(t_fine, x_exact_fine(1,:), 'k-', 'LineWidth', 3, 'DisplayName', 'Exact Solution');

% Plot the numerical solutions
for i = 1:length(results)
    plot(results(i).t, results(i).x(1,:), '--', 'LineWidth', 1.5, ...
        'DisplayName', ['dt = ', num2str(results(i).dt)]);
end
hold off;
title('RK4 Grid Independence Study: Lateral Speed');
xlabel('Time (s)');
ylabel('Yaw Rate (rad/s)');
legend('show', 'Location', 'southeast');
grid on;

figure;
hold on;
% Plot the exact solution (Yaw Rate)
t_fine = 0:0.001:T;
x_exact_fine = f_exact(t_fine);
plot(t_fine, x_exact_fine(1,:), 'k-', 'LineWidth', 3, 'DisplayName', 'Exact Solution');

% Plot the numerical solutions
for i = 1:length(results)
    plot(results(i).t, results(i).x(1,:), '--', 'LineWidth', 1.5, ...
        'DisplayName', ['dt = ', num2str(results(i).dt)]);
end
hold off;
title('RK4 Grid Independence Study: Yaw Rate');
xlabel('Time (s)');
ylabel('Yaw Rate (rad/s)');
legend('show', 'Location', 'southeast');
grid on;



% 4. Create a log-log plot for visual verification of convergence
figure;
plot(log(step_sizes), log(errors), 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Numerical Error');
hold on;
% Plot a reference line with a slope of 4
C = errors(1) / (step_sizes(1)^4); % Estimate constant C for the reference line
ref_line = C * step_sizes.^4;
plot(log(step_sizes), log(ref_line), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Reference Slope of 4');
hold off;

title('RK4 Convergence: Error vs. Step Size');
xlabel('Step Size (h)');
ylabel('True Error at T=5s');
legend('show', 'Location', 'northwest');
grid on;
axis tight;


% ==== FUNCTION DEFINITIONS ====

% Generic IVP solver
function [t, y] = solveIVP(f, tspan, y0, h, solver)
    % Ensure tspan(2) is a multiple of h for accurate endpoint
    num_steps = round((tspan(2) - tspan(1)) / h);
    t = tspan(1):h:(tspan(1) + num_steps * h);
    % Make sure final time point is included, corrects for floating point issues
    if abs(t(end) - tspan(2)) > 1e-9 
        t = [t, tspan(2)]; 
    end
    y = zeros(length(y0), length(t));
    y(:,1) = y0;
    for n = 1:length(t) - 1
        h_step = t(n+1) - t(n); % Use actual step in case last one is smaller
        y(:,n+1) = solver(f, t(n), y(:,n), h_step);
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
