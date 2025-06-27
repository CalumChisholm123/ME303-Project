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
        L2_norm_lat = sqrt(sum(actLatSpeed - latSpeedValues).^2);
        L2_norm_yaw = sqrt(sum(actYawRates - yawRateValues).^2);
        
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
