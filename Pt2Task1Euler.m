%Euler Method (forward)
function [t,x] = euler_1(ode, tspan, x0, h)

    t = tspan(1) : h : tspan(2);
    n = length(t);
    x = zeros(length(x0), n);
    x(:,1) = x0;

    for i = 1:length(t) -1;
        x(:, i+1) = x(:, i) + h * ode(t(i), x(:, i));
    end
end

stepSize = [0.1, 0.05, 0.01, 0.001];

figure;
hold on;
for h = stepSize
    [t, x] = euler_1(bicycle_ODE, tspan, x0, h);
    plot(t, x(2,:), 'DisplayName', ['h = ', num2str(h)]);
end

title('Euler Method - Grid Independence Check');
xlabel('Time (s)');
ylabel('Yaw rate (rad/s)');
legend;
grid on;

h_selected = 0.01; % Chosen based on grid study
[t, x] = euler_method(bicycle_ode, tspan, x0, h_selected);

figure;
plot(t, x(2,:), 'b', 'LineWidth', 1.5);
title('Final Euler Solution (h = 0.01)');
xlabel('Time (s)');
ylabel('Yaw rate (rad/s)');
grid on;