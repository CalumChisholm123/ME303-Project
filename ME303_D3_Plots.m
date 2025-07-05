% RK4 method solver for vehicle dynamics - Grid independence study

clear;
clc; 
close all;

addpath('/')

% Parameters
m = 1400;          % kg
a = 1.14;          % m
    b = 1.33;          % m
Caf = 25000;       % N/rad
Car = 21000;       % N/rad
Iz = 2420;         % kg·m^2
u = 75 / 3.6;      % Convert from km/h to m/s
delta = 0.1;       % Step steering input (rad)

m_list = 650:150:m;


% Time setup
T = 5;             % Total simulation time (s)
dt = [0.1, 0.05, 0.02, 0.01, 0.005]; % Step Values
h = 1e-5

x0 = [0,0];

figure;
hold on; 
colors = lines(length(m_list)); % Generate distinct colors for each line
legend_entries = {}; % Initialize a cell array for legend entries

% --- Main Simulation Loop ---
% Iterate through each mass value in m_list
for n = 1:length(m_list)
    m = m_list(n); % Current mass for this iteration
    
    % Define the A matrix (System Matrix) for the current mass
    A = [- (Caf + Car)/(m*u),   (-a*Caf + b*Car)/(m*u) - u;
         (-a*Caf + b*Car)/(Iz*u), -(a^2*Caf + b^2*Car)/(Iz*u)];
    
    B = [Caf/m;
         a*Caf/Iz];
    

    f = @(t, x) A * x + B * delta;
    

    [t, x] = solveIVP(f, [0,5], x0, h, @rk4);
    
    plot(t, x(1,:), 'LineWidth', 2, 'Color', colors(n,:));
    
    % --- Legend ---
    % Create a descriptive entry for the legend
    legend_entries{n} = ['Mass = ', num2str(m), ' kg'];
end

% --- Finalize the Plot ---
hold off; % Release the plot hold
grid on;  % Add a grid for better readability

% Add labels and a title
title('Effect of Mass on Vehicle Lateral Speed', 'FontSize', 14);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Lateral Speed (m/s)', 'FontSize', 12);

% Add the legend to the plot
legend(legend_entries, 'Location', 'SouthEast');


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
        L2_norm_lat = sqrt(sum((actLatSpeed - latSpeedValues).^2));
        L2_norm_yaw = sqrt(sum((actYawRates - yawRateValues).^2));

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
