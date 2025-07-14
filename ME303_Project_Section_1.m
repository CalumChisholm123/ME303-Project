%% ==================== 2.3 Euler Method ====================

syms x
f = exp(x);

% Calculate the Taylor series expansions
taylor1 = taylor(f, 'Order', 2);  
taylor3 = taylor(f, 'Order', 4);
taylor5 = taylor(f, 'Order', 6);
taylor10 = taylor(f, 'Order', 11);

figure(1);
hold on;

fplot(f, [0, 5], 'LineWidth', 2, 'DisplayName', 'e^x');

fplot(taylor1, [0, 5], '--', 'DisplayName', 'Order 1');
fplot(taylor3, [0, 5], '--', 'DisplayName', 'Order 3');
fplot(taylor5, [0, 5], '--', 'DisplayName', 'Order 5');
fplot(taylor10, [0, 5], '--', 'DisplayName', 'Order 10');

grid on;
title('Taylor Series Approximations of e^x');
xlabel('X Values');
ylabel('Y Values');
xlim([0, 5]);
ylim([0, 150]); 
legend show; 
hold off; 


clear;
clc;
close all;

% Calculate the Error 
clear;
clc;
close all;

syms x
f = exp(x);

figure;
hold on;

orders_to_plot = [1, 3, 5, 10]; 

for i = 1:length(orders_to_plot)
    n = orders_to_plot(i);
    

    taylor_approx = taylor(f, 'Order', n + 1);
    

    error_func = abs(f - taylor_approx);
    

    fplot(error_func, [0, 5], 'DisplayName', ['Order ' num2str(n)]);
end


grid on;
title('Absolute Error for Series Approximations of e^x');
xlabel('x');
ylabel('Absolute Error');
legend show;
hold off;

%% ==================== 2.4 Eulers Method Grid Spacing ====================
clear; clc;

f = @(x, y) exp(x); 

x0 = 0; 
y0 = 0; 

x_end = 5;  

dx_values = [0.1, 0.05, 0.005]; 

figure;
hold on; 

legend_entries = {};

for i = 1:length(dx_values)
    [x, y] = euler_method(f, [x0, x_end], y0, dx_values(i));
    plot(x, y, 'LineWidth', 1.5); 
    legend_entries{i} = ['Euler, dx = ', num2str(dx_values(i))];
end

x_exact = linspace(x0, x_end, 1000);
y_exact = exp(x_exact) - 1;  
plot(x_exact, y_exact, 'k:', 'LineWidth', 2);
legend_entries{end+1} = 'Exact solution';


hold off;
xlabel('x');
ylabel('y(x)');
title('Solution of dy/dx = e^x using Euler''s Method');
grid on;
legend(legend_entries, 'Location', 'northwest');
xlim([x0 x_end]);
ylim([0 20]); 


function [x, y] = euler_method(f, xspan, y0, dx)
    x = xspan(1):dx:xspan(2);  
    y = zeros(size(x));         
    y(1) = y0;                
    
    for i = 1:length(x)-1
        y(i+1) = y(i) + dx * f(x(i), y(i));
    end
end
