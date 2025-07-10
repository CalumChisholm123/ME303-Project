classdef VehicleDynamics
    % VehicleDynamics Class for simulating vehicle dynamics using RK4
    
    properties
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
        dt = 0.01;         % Step size
        tspan = [0, 30];   % Time span
    end
    
    methods
        % Method to define A and B matrices
        function [A, B] = getMatrices(obj)
            A = [- (obj.Caf + obj.Car)/(obj.m * obj.u), (-obj.a * obj.Caf + obj.b * obj.Car)/(obj.m * obj.u) - obj.u;
                 (-obj.a * obj.Caf + obj.b * obj.Car)/(obj.Iz * obj.u), - (obj.a^2 * obj.Caf + obj.b^2 * obj.Car)/(obj.Iz * obj.u)];
            B = [obj.Caf / obj.m; obj.a * obj.Caf / obj.Iz];
        end
        
        % Method to define ODE function
        function f = getODEFunction(obj, A, B)
            f = @(t, x) A * x + B * obj.delta;
        end
        
        % Generic IVP solver
        function [t, y] = solveIVP(~, f, tspan, y0, h, solver)
            t = tspan(1):h:tspan(2);
            y = zeros(length(y0), length(t));
            y(:, 1) = y0;
            for n = 1:length(t) - 1
                y(:, n + 1) = solver(f, t(n), y(:, n), h);
            end
        end
        
        % Runge-Kutta 4th order method
        function ynew = rk4(~, f, t, y, h)
            k1 = f(t, y);
            k2 = f(t + 0.5 * h, y + 0.5 * h * k1);
            k3 = f(t + 0.5 * h, y + 0.5 * h * k2);
            k4 = f(t + h, y + h * k3);
            ynew = y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
        end
        
        % Method to run the simulation
        function [t, x] = runSimulation(obj)
            % Get A and B matrices
            [A, B] = obj.getMatrices();
            
            % Define ODE function
            f = obj.getODEFunction(A, B);
            
            % Initial condition
            x0 = [0; 0];
            
            % Solve IVP using RK4
            [t, x] = obj.solveIVP(f, [0, obj.T], x0, obj.dt, @obj.rk4);
        end
    end
end