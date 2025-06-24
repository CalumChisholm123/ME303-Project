%Parameters for Vehicle
m = 1400;          % kg
a = 1.14;          % m
b = 1.33;          % m
Caf = 25000;       % N/rad
Car = 21000;       % N/rad
Iz = 2420;         % kg·m^2
u = 75 / 3.6;      % Convert from km/h to m/s
delta = 0.1;       % Step steering input (rad)

%Matrix
A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
     (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];

B = [Caf*delta/m, a*Caf*delta/Iz];

A_u = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
       (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];

B_u = [Caf*delta/m, a*Caf*delta/Iz];

%Define IVP
tspan = [0,5];
y0 = 1;
x0 = [0 0];
stepSizes = [0.1, 0.05, 0.01, 0.001];

h_selected = 0.01; % Chosen based on grid study

%Plot ydots
% figure(1)
% hold on
% for i = 1 : length(stepSizes)
%     h = stepSizes(i);
%     [t, y_dot, yawrate] = fwd_euler(A, B, tspan, h);
%     plot(t, y_dot, 'DisplayName', ['h = ', num2str(h)])
%     hold on
% end
% title('Euler Method - Lateral velocity vs Time');
% xlabel('Time (s)');
% ylabel('$\dot{y}$ (km/h)', 'Interpreter', 'latex');
% legend;
% grid on;
% 
% figure(2)
% hold on
% for i = 1 : length(stepSizes)
%     h = stepSizes(i);
%     [t, y_dot, yawrate] = fwd_euler(A, B, tspan, h);
%     plot(t, yawrate, 'DisplayName', ['h = ', num2str(h)])
%     hold on
% end
% title('Euler Method - Yaw Rate - Grid Independence Check');
% xlabel('Time (s)');
% ylabel('Yaw rate (rad/s)');
% legend;
% grid on;

function [t,x,y] = fwd_euler(A, B, tspan, h)
    t = tspan(1) : h : tspan(2);
    n = length(t);
    x = zeros(1,n); %ydot
    y = zeros(1,n); %phidot
    x(1) = 0;
    y(1) = 0;
    for i = 1 : length(t) - 1
        x(i+1) = x(i) + (A(1,1)*x(i) + A(1,2) * y(i) + B(1))*h;
        y(i+1) = y(i) + (A(2,1)*x(i) + A(2,2) * y(i) + B(2))*h;
    end
end

function [t, x, y] = B1(tspan, h, u)
    % u is longitudinal speed input in km/h -> convert to m/s in the func.
    u = u ./ 3.6; %convert km/h to m/s;
    
    
end
