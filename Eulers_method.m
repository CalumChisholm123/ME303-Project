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

% Grid convergence
times = tspan(1):1:tspan(2);
actLatSpeed = -13.0964.*exp(-1.9745.*times) + 24.4684*exp(-0.9839.*times) - 11.3720;
actYawRates = -0.2496.*exp(-1.9745.*times) - 0.6962.*exp(-0.9839.*times) + 0.9457;

[L2_lat_speed, L2_yaw_rate] = grid_check_L2_norm(f, tspan, x0, stepSizes, actLatSpeed, actYawRates);


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

function [t,x] = euler_1(ode, tspan, x0, dt)

    t = tspan(1) : dt : tspan(2);
    n = length(t);
    x = zeros(length(x0), n);
    x(:,1) = x0;

    for i = 1:length(t) -1;
        x(:, i+1) = x(:, i) + dt .* ode(t(i), x(:, i));
    end
end

function [L2LatSpeedNorm, L2YawRateNorm] = grid_check_L2_norm(ODE, tspan, x0, stepSizes, actLatSpeed, actYawRates)
    L2LatSpeedNorm = zeros(length(stepSizes)); %Store L2 values for different stepSizes
    L2YawRateNorm = zeros(length(stepSizes));
    times = tspan(1):1:tspan(2);

    for i = 1:length(stepSizes)
        [t, x_dot] = euler_1(ODE, tspan, x0, stepSizes(i));
        latSpeedValues = zeros(0, length(times)); %vector to store length(times) values
        yawRateValues = zeros(0, length(times));

        for j = 1:length(times)
            temp_indx = find(t == times(j));
            latSpeedValues(j) = x_dot(1,temp_indx);
            yawRateValues(j) = x_dot(2, temp_indx);
        end 
        L2_norm_lat = sqrt(sum(actLatSpeed - latSpeedValues).^2);
        L2_norm_yaw = sqrt(sum(actYawRates - yawRateValues).^2);
        
        L2LatSpeedNorm(i) = L2_norm_lat;
        L2YawRateNorm(i) = L2_norm_yaw;
    end
    
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

