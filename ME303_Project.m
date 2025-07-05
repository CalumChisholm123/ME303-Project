% RK4 method solver for vehicle dynamics - Grid independence study

clear all;
clc

% Parameters
m = 1400;          % kg
a = 1.1859;          % distance of CM to front axle 
b = 1.28;          % distance of CM to rear axle 
Caf = 25000;       % N/rad
Car = 21000;       % N/rad
Iz = 2420;         % kg·m^2
delta = 0.1;       % Step steering input (rad)

% Time setup
tspan = [0,5];
T = 5;             % Total simulation time (s)
dt = 0.001;% Step Values 

N = T/dt; % Number of Steps 
t = linspace(0,T,N+1); % Time Vector 

% Define A and B matrices (constant for constant u)
%A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
%     (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];

B = [Caf/m; a*Caf/Iz];

% Define ODE function: dx/dt = A*x + B*delta
f = @(t, x) A * x + B * delta;

% Initial condition: lateral velocity and yaw rate both zero
x0 = [0; 0];

%[t, x] = solveIVP(f, [0, T], x0, dt, @rk4);

% figure(1)
% plot(t,x(2,:),"b", LineWidth=2);
% hold on;
% title("Time vs Yaw - RK4")
% ylabel("Yaw Rate (degrees)");
% xlabel("Time");
% grid on;
% 
% figure(2)
% plot(t,x(1,:),"b", LineWidth=2);
% hold on;
% title("Time vs Lateral Velocity - RK4")
% ylabel("Lateral Velocity (m/s)");
% xlabel("Time");
% grid on;
% 
% 

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

%Euler Method (forward)
function [t,x] = euler_1(ode, tspan, x0, dt)

    t = tspan(1) : dt : tspan(2);
    n = length(t);
    x = zeros(length(x0), n);
    x(:,1) = x0;

    for i = 1:length(t) -1
        x(:, i+1) = x(:, i) + dt * ode(t(i), x(:, i));
    end
end


% u_values = (1:0.0001:5000);
% 
% for i = 1:length(u_values)
%     u = u_values(i);
% 
%     A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
%      (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];
% 
% 
%     lambda = eig(A);
%     disp(lambda)
% 
%     if any(real(lambda) >= 0)
%         u = u*3.6;
%         disp(u)
% 
%         disp('UNSTABLE: At least one eigenvalue has a non-negative real part.');
% 
%         break;
%     end
% 
% end

stability_factor = Car * b - Caf * a;

if stability_factor == 0
    disp('The bicycle is neutrally steered; critical velocity is theoretically infinite.');
    u_crit = inf;
else
    % Calculate the term inside the square root
    radicand = (-Caf * Car * (a + b)^2) / (m * stability_factor);

    % Check if the bicycle is understeering or oversteering
    if radicand >= 0
        % Oversteering case: Calculate critical velocity
        u_crit = sqrt(radicand);
        fprintf('The critical velocity is: %.2f m/s\n', u_crit);
    else
        % Understeering case: No real critical velocity
        disp('The bicycle is understeering and stable at all speeds. Critical velocity is not real.');
        u_crit = NaN; % Assign Not-a-Number for non-real results
    end
end



% Calums Speed Changer 
% u_values = [20,50,75,100,200,300] / 3.6;
% 
% figure(1);
% hold on;
% for i = 1:length(u_values)
%     u = u_values(i);
% 
%     A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
%      (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];
% 
%     B = [Caf/m; a*Caf/Iz];
% 
%     f = @(t, x) A * x + B * delta;
%     [t, x] = solveIVP(f, [0, T], x0, dt, @rk4);
% 
%     plot(t,x(1,:),'Displayname', ['u = ',num2str(u)]);
% 
% end  
% legend;
% xlabel('Times (s)');
% ylabel('Lateral Velocity (m/s)')
% title("Lateral Velocity Different u values - RK4");
% clear all;

% figure(2);
% hold on;
% for p = 1:length(u_values)
%     u = u_values(p);
% 
%     A = [- (Caf + Car)/(m*u), (-a*Caf + b*Car)/(m*u) - u;
%      (-a*Caf + b*Car)/(Iz*u), - (a^2*Caf + b^2*Car)/(Iz*u)];
% 
%     B = [Caf/m; a*Caf/Iz];
% 
%     f = @(t, x) A * x + B * delta;
% 
%     [t, x] = euler_1(f, tspan, x0, dt);
% 
%     plot(t,x(1,:),'Displayname', ['u = ',num2str(u)]);
% 
% end  
% legend;
% title("Different u values - Eulers");