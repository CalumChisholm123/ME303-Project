clear all;
clc

% Parameters
m = 1400;          % kg
a = 1.14;          % m
b = 1.33;          % m
Caf = 25000;       % N/rad
% Car = 21000;       % N/rad
Iz = 2420;         % kg·m^2
u = 75 / 3.6;      % Convert from km/h to m/s
delta = 0.1;       % Step steering input (rad)
x0 = [0; 0];
Car_span = [10000, 50000];

%[C_rears, stability_list, lambdas_real1, lambdas_real2, lambdas_imag1, lambdas_imag2] = check_stability(Caf, Car_span, 1, m, u, a, b, Iz);

%Plotting the real and imaginary components and applying a gradient
max_C_rear = Car_span(2);
% plot_gradient(C_rears, lambdas_real1, lambdas_real2, lambdas_imag1, lambdas_imag2)
%Plot real lambdas
% figure(1)
% plot(C_rears, lambdas_real1, "DisplayName", "\lambda_1")
% hold on
% plot(C_rears, lambdas_real2, "DisplayName", "\lambda_2")
% title("Real \lambda_1 and \lambda_2 for constant C_{\alphaf} = 25000 N/rad and variable C_{\alphar}")
% xlabel("C_{\alphar}")
% ylabel("\lambda")
% legend
% 
% figure(2)
% plot(C_rears, lambdas_imag1, "DisplayName", "\lambda_1")
% hold on
% plot(C_rears, lambdas_imag2, "DisplayName", "\lambda_2")
% title("Imaginary \lambda_1 and \lambda_2 components for constant C_{\alphaf} = 25000 N/rad and variable C_{\alphar}")
% xlabel("C_{\alphar}")
% ylabel("\lambda imaginary")
% legend

Cars_list = [18000 21000 21200 30000];
tspan = [0, 5];
[times, lat_vel_results, yaw_rate_results] = solve_diff_Cars(Cars_list, tspan, 0.001, x0, m, u, a, b, Iz, Caf, delta);

% Function
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

%Check all C_alpha values, also solve the ODE for each C_alpha? No need...
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

%THIS FUNCTION DID NOT END UP GETTING USED.
function plot_gradient(C_rears, lambdas_real1, lambdas_real2, lambdas_imag1, lambdas_imag2)
    n = length(lambdas_real1);
    cmap = winter(n);
    figure;
    hold on;
    for i = 1:n
        plot(lambdas_real1(i), lambdas_imag1(i), '.', 'Color', cmap(i,:));
        hold on
        plot(lambdas_real2(i), lambdas_imag2(i), '.', 'Color', cmap(i,:));
        title("Eigven values for constant C_{\alphaf} = 25000 N/rad and variable C_{\alphar}")
        xlabel("Real component")
        ylabel("Imaginary component")
        grid on
    end
end

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