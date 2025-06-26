%% Project Skeleton Code

clear; clc; close all;

%%initial parameter: unit: m, degree, rad/sec
% Find r4 and r6
r1 = 5.0*10^(-2); % m  o2o3
r2 = 2.0*10^(-2); % m  o2a2
r3 = 8.5*10^(-2); % m  o3b
r5 = 3.0*10^(-2); % m  cb
r7 = 10.0*10^(-2); % m  o3o4

theta2 = 0:1:360; % from 0 to 360 with step 1: [0,1,2,3,4....360]
dtheta2 = 10;
ddtheta2 = 0; 

% TIPS:  

% cosd(x) - is a cosine of x, where x in degrees
% cos(x) - is a cosine of x, where x in radians
% using '.*' enables element-wise multiplication
% accordingly, '.^' element-wise exponent
% [a1 a2 a3].^[b1 b2 b3] = [a1*b1 a2*b2 a3*b3]
% '*' is matrix multiplication

%% Part 1- Calculations for kinematic variables, caculated based on loop closure eqn

% Hint: Check this for all other angles too

% r4, theta3 caculated based on loop eqn 1
r4 = sqrt(r1.^2 + r2.^2 - 2.*r2.*r1.*cosd(theta2));

% Hint: Check if the angle needs to be adjusted to its true value
theta3 = atand(r2 .* sind(theta2)./ (r1 - r2 .* cosd(theta2))); 
theta3 = 180-theta3; 

% theta5, r6 caculated based on loop eqn 2 
% Hint: Check if the angle needs to be adjusted to its true value
theta5 = acosd((-r7-r3.*cosd(theta3))./r5);

r6 = r3.*sind(theta3) + r5.*sind(theta5);

%% Take time derivative of loop eqn (d/dt) 
% and solve them for dtheta3, dtheta5 & dr6
% and the same for the second derivatives. 

dr4 = r2 .* dtheta2 .* sind(theta2 - theta3);

dtheta3 = (r2 .* dtheta2 ./ r4) .* cosd(theta2-theta3);

dtheta5 = (-r3 .* dtheta3 .* sind(theta3)) ./ (r5 .* sind(theta5));

dr6 = r3 .* dtheta3 .* cosd(theta3) + r5 .* dtheta5 .* cosd(theta5);

ddtheta3 = (r2 ./ r4) .* dtheta2 .* (dtheta3 - dtheta2) .* sind(theta2 - theta3) - (r2 .* dr4 ./ r4) .* dtheta2 .* cosd(theta2 - theta3); 

ddtheta5 = (-r3 .* dtheta3.^2 .* cosd(theta3) -r3 .* ddtheta3 .* sind(theta3) -r5 .* dtheta5.^2 .* cosd(theta5)) ./ (r5 .* sind(theta5));

ddr6 = r3 .* (ddtheta3 .* cosd(theta3) - (dtheta3.^2).*sind(theta3)) + r5 .* (ddtheta5 .* cosd(theta5) - (dtheta5.^2).*sind(theta5));

ddr3 = 0;
%% Plot vars;

% Plot all desired deliverables. 
figure (1)
plot(theta2,theta3)
grid on;
title('$\theta_3$ vs $\theta_2$', 'Interpreter','latex')
xlabel('\theta_2   unit: degree')
ylabel('\theta_3   unit: degree')

figure (2)
plot(theta2,theta5)
grid on;
title('$\theta_5$ vs $\theta_2$', 'Interpreter','latex')
xlabel('\theta_2   unit: degree')
ylabel('\theta_5   unit: degree')

figure (3)
plot(theta2, r6)
grid on;
title('$r_6$ vs $\theta_2$', 'Interpreter','latex')
xlabel('\theta_2   unit: degree')
ylabel('r_6  unit: cm')

figure (4)
plot(theta2, r4)
grid on;
title('$r_4$ vs $\theta_2$', 'Interpreter','latex')
xlabel('\theta_2   unit: degree')
ylabel('r_4  unit: cm')
 
figure (5)
plot(theta2,dtheta3)
grid on;
title('$\dot{\theta_3}$ vs $\theta_2$', 'Interpreter','latex')
xlabel('\theta_2   unit: degree')
ylabel('$\dot{\theta_3}$  unit: degree/sec', 'Interpreter','latex')

figure (6)
plot(theta2,dtheta5)
grid on;
title('$\dot{\theta_5}$ vs $\theta_2$', 'Interpreter','latex')
xlabel('\theta_2   unit: degree')
ylabel('$\dot{\theta_5}$  unit: degree/sec', 'Interpreter','latex')

figure (7)
plot(theta2, dr4)
grid on;
title('$\dot{r_4}$ vs $\theta_2$', 'Interpreter','latex')
xlabel('\theta_2   unit: degree')
ylabel('$\dot{r_4}$  unit: cm/sec', 'Interpreter','latex')

figure (8)
plot(theta2, dr6)
grid on;
title('$\dot{r_6}$ vs $\theta_2$', 'Interpreter','latex')
xlabel('\theta_2   unit: degree')
ylabel('$\dot{r_6}$  unit: cm/sec', 'Interpreter','latex')


figure (9)
plot(theta2, ddtheta3)
grid on;
title('$\ddot{\theta_3}$ vs $\theta_2$', 'Interpreter','latex')
xlabel('\theta_2   unit: degree')
ylabel('$\ddot{\theta_3}$  unit: degree/sec$^2$', 'Interpreter','latex')


figure (10)
plot(theta2, ddtheta5)
grid on;
title('$\ddot{\theta_5}$ vs $\theta_2$', 'Interpreter','latex')
xlabel('\theta_2   unit: degree')
ylabel('$\ddot{\theta_5}$  unit: degree/sec$^2$', 'Interpreter','latex')


figure (11)
plot(theta2, ddr6)
grid on;
title('$\ddot{r_6}$ vs $\theta_2$', 'Interpreter','latex')
xlabel('\theta_2   unit: degree')
ylabel('$\ddot{r_6}$  unit: cm/sec$^2$', 'Interpreter','latex')
% *****************************************************

%% Part 2 - Force and Moment Calculation

%initial parameter: unit: m, degree, rad/sec
r1 = 5.0*10^(-2); % m  o2o3
r2 = 2.0*10^(-2); % m  o2a2
r3 = 8.5*10^(-2); % m  o3b
r5 = 3.0*10^(-2); % m  cb
r7 = 10.0*10^(-2); % m  o3o4

dtheta2 = -10; % theta2 dot
ddtheta2 = 0; % theta2 doble-dot - second derivative

rho = 2700; % kg/m3
d = 0.5*10^(-2); % m

m2 = 3.14 * 0.25 * d^2 * r2 * rho ; % link 2, o2a2 kg
m3 = 3.14 * 0.25 * d^2 * r3 * rho; % link 3, o3b, kg
m5 = 3.14 * 0.25 * d^2 * r5 * rho; % link 5, bc, kg
m9 = 5 * 10^(-3); % slider 6, kg 
m10 = 5 * 10^(-3); % slider 4, kg 

I_G3 = m3*r3^2/12;
I_G5 = m5*r5^2/12;
%%
M12_list = [];
theta2_list = [];
Fs_list = [];  % shaking force
alpha_s_list = []; % direction of a shaking force

% Shaking moment
Ms_list =[]; 

% Forces:
% Ground, link2, slider 4 
F12_list = [];
F102_list = [];

% Ground, link3, slider 4 
F13_list = [];
F53_list = [];
F103_list = [];

% Ground, link5, slider 6 
F59_list = [];
F19_list = [];

% Angles at which forces are acting:
F12_alpha = [];
F102_alpha = [];
F13_alpha = [];
F53_alpha = [];
F59_alpha = [];

% Indeces of the points to be marked up:
marker_indices = [];

for theta2 = 0:1:360

    % kinematic variables are calculated based on loop eqn
    r4 = sqrt(r1.^2 + r2.^2 - 2.*r2.*r1.*cosd(theta2));
    
    theta3 = atand(r2 .* sind(theta2)./ (r1 - r2 .* cosd(theta2))); 
    theta3 = 180-theta3; 

    theta5 = acosd((-r7-r3.*cosd(theta3))./r5);

    r6 = r3.*sind(theta3) + r5.*sind(theta5);
    
    dtheta3 = (r2 .* dtheta2 ./ r4) .* cosd(theta2-theta3);
    
    dtheta5 = (-r3 .* dtheta3 .* sind(theta3)) ./ (r5 .* sind(theta5));
    
    dr4 = r2 .* dtheta2 .* sind(theta2 - theta3);
    
    ddtheta3 = (r2 ./ r4) .* dtheta2 .* (dtheta3 - dtheta2) .* sind(theta2 - theta3) - (r2 .* dr4 ./ r4) .* dtheta2 .* cosd(theta2 - theta3); 
     
    ddtheta5 = (-r3 .* dtheta3.^2 .* cosd(theta3) -r3 .* ddtheta3 .* sind(theta3) -r5 .* dtheta5.^2 .* cosd(theta5)) ./ (r5 .* sind(theta5));

    dr6 = r3 .* dtheta3 .* cosd(theta3) + r5 .* dtheta5 .* cosd(theta5);

    ddr6 = r3 .* (ddtheta3 .* cosd(theta3) - (dtheta3.^2).*sind(theta3)) + r5 .* (ddtheta5 .* cosd(theta5) - (dtheta5.^2).*sind(theta5));


    % Get force vector
    B = get_ma_vector(m2, m3, m5, m9, m10, ...
        r2, r3, r5, ...
        theta2, theta3, theta5, ...
        dtheta2, dtheta3, dtheta5, ...
        ddtheta3, ddtheta5, ...
        ddr6, ...
        I_G3, I_G5);
    
    % Get acceleration matrix 
    A = get_A_matrix(...
        r2,r3,r4,r5,...
        theta2,theta3,theta5);

    x = A\ B; % Ax = B, solution for x; note that in MATLAB: A\B = B/A
    
    % M12:
    M12 = x(11);
    M12_list = [M12_list; M12];
    
    F12x = x(1);
    F12y = x(2);
    F13x = x(5);
    F13y = x(6);
    
    Fs = sqrt((F12x + F13x)^2 + (F12y + F13y)^2); % Shaking force
    Fs_list = [Fs_list; Fs];
    
    
    % Directions of all forces:    
    fx = (F12x + F13x);
    fy = (F12y + F13y);
    alpha_s = atan(fx\fy);
    if fx < 0
         alpha_s = alpha_s + pi;
    end 
    alpha_s_list = [alpha_s_list,alpha_s];

    Ms = - M12 +  r1 .* F13y; % Shaking moment
    Ms_list = [Ms_list; Ms];
    
    
    % Magnitudes of all forces: 
    F12_list = [F12_list; sqrt(x(1)^2 + x(2)^2)];
    F102_list = [F102_list; sqrt(x(3)^2 + x(4)^2)];
    F13_list = [F13_list; sqrt(x(5)^2 + x(6)^2)];
    F53_list = [F53_list; sqrt(x(7)^2 + x(8)^2)];
    F59_list = [F59_list; sqrt(x(9)^2 + x(10)^2)];
    F103_list = [F103_list; sqrt(x(12)^2)];
    F19_list = [F19_list; sqrt(x(13)^2)];
    
    
    % Directions of all forces:    
    fx = x(1);
    fy = x(2);
    alpha_f = atan(fx\fy);
    if fx < 0
        alpha_f = alpha_f + pi;
    end
    F12_alpha = [F12_alpha; alpha_f];

    fx = x(3);
    fy = x(4);
    alpha_f = atan(fx\fy);
    if fx < 0
        alpha_f = alpha_f + pi;
    end
    F102_alpha = [F102_alpha; alpha_f];
         
    fx = x(5);
    fy = x(6);
    alpha_f = atan(fx\fy);
    if fx < 0
        alpha_f = alpha_f + pi;
    end
    F13_alpha = [F13_alpha; alpha_f];
    
    fx = x(7);
    fy = x(8);
    alpha_f = atan(fx\fy);
    if fx < 0
        alpha_f = alpha_f + pi;
    end
    F53_alpha = [F53_alpha; alpha_f];
    
  
    fx = x(9);
    fy = x(10);
    alpha_f = atan(fx\fy);
    if fx < 0
        alpha_f = alpha_f + pi;
    end
    F59_alpha = [F59_alpha; alpha_f];
    
    % Collecting the values of theta2:
    theta2_list = [theta2_list, theta2];
     
    % Collect the indices to be marked up in polar plots:
    if rem(theta2,30) == 0 % for each 30 degrees
        marker_indices = [marker_indices; theta2];
    end
   
    
end


% Plots:


figure (13)
plot(theta2_list,M12_list)
grid on;
title('M_{12} vs \theta_2')
xlabel('\theta_2   unit: degree')
ylabel('M12   unit: N-m')




figure (26)
plot(theta2_rad,Ms_list)
grid on;
title('M_s vs \theta_2')
xlabel('\theta_2   unit: radian')
ylabel('Ms   unit: N-m')


