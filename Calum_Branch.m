% % Generic IVP solver
% function [t, y] = solveIVP(f, tspan, y0, h, solver)
%     t = tspan(1):h:tspan(2);
%     y = zeros(length(y0), length(t));
%     y(:,1) = y0;
%     for n = 1:length(t) - 1
%         y(:,n+1) = solver(f, t(n), y(:,n), h);
%     end
% end
% 
% % Runge-Kutta 4th order method
% function ynew = rk4(f, t, y, h)
%     k1 = f(t, y);
%     k2 = f(t + 0.5 * h, y + 0.5 * h * k1);
%     k3 = f(t + 0.5 * h, y + 0.5 * h * k2);
%     k4 = f(t + h, y + h * k3);
%     ynew = y + (h / 6) * (k1 + 2*k2 + 2*k3 + k4);
% end
% 
% 
% clear all;
% clc
% 
% syms x
% f = exp(x);
% 
% % Calculate the Taylor series expansions
% taylor1 = taylor(f, 'Order', 2);  % Order n means the polynomial has degree n-1
% taylor3 = taylor(f, 'Order', 4);
% taylor5 = taylor(f, 'Order', 6);
% taylor10 = taylor(f, 'Order', 11);
% 
% % Set up the plot
% figure(1); % Create a new figure window
% hold on; % Hold the plot to overlay multiple graphs
% 
% % Plot the original function using fplot
% fplot(f, [0, 5], 'LineWidth', 2, 'DisplayName', 'e^x');
% 
% % Plot the Taylor series approximations using fplot
% fplot(taylor1, [0, 5], '--', 'DisplayName', 'Order 1');
% fplot(taylor3, [0, 5], '--', 'DisplayName', 'Order 3');
% fplot(taylor5, [0, 5], '--', 'DisplayName', 'Order 5');
% fplot(taylor10, [0, 5], '--', 'DisplayName', 'Order 10');
% 
% % Add plot details
% grid on;
% title('Taylor Series Approximations of e^x');
% xlabel('X Values');
% ylabel('Y Values');
% xlim([0, 5]);
% ylim([0, 150]); % Adjust y-axis for better visualization
% legend show; % Display the legend
% hold off; % Release the plot hold

% 
clear all;
clc
f = @(x,y) exp(x);  % Also fixed function definition syntax
x0 = 0;
y0 = 1;
xfinal = 5;

h = 0.05;

N = (xfinal -x0)/h;

x(1) = x0;
y(1) = y0;

for n = 1:N
    y(n+1) = y(n) + h*f(x(n),y(n));
    x(n+1) = x(n) + h;
end

figure(1)
grid on;
plot(x,y)

title('Numerical Solution of IVP using Euler Method, Step Size 0.05');
xlabel("X Values")
ylabel("Y Values")    