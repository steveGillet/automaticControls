b = 1;
m1 = 10; % [kg]
m2 = 350; % [kg]
kw = 500000; % [Nm]
ks = 10000; % [Nm]

% % Original Suspension
% A = [-b/m2 -1/m2 b/m2 0; ks 0 -ks 0; b/m1 1/m1 -b/m1 -1/m1; 0 0 kw 0]; 
% B = [0; 0; 0; -kw];       
% C = [1 0 0 0];        
% D = [0];  

% Automatic Suspension
A = [0 0 0; 0 0 -1/m1; 0 kw 0]; 
B = [0 -1/m2; 0 1/m1; -kw 0];       
C = [1 0 0];        
D = [0];  

carSusSS = ss(A, B, C, D);

carSusTF = tf(carSusSS);
disp(carSusTF);

% plot(wheelDisplacement.time, wheelDisplacement.Data, 'b-', ...
%      carDisplacement.time, carDisplacement.Data, 'r-');

% xlabel('Time [s]');
% ylabel('Displacement [m]');
% legend('Wheel Displacement', 'Car Displacement');
% title(sprintf('b = %d, Wheel vs Car Displacement',b));
% grid on;

% Auto Sus Simulation

t = 0:0.01:10; 
r = zeros(size(t)); 
r(t >= 1) = 1;
% f = zeros(size(t)); 
% f(t >= 1) = 100;
f = 1 * sin(35.6 * 2 * pi * t);
u = [r; f];

[y, t, x] = lsim(carSusSS, u, t);
format longG
disp(eig(A));

% Plot the states
figure;
plot(t, x(:,1), 'b-', t, x(:,2), 'r-', t, r, 'g-');
xlabel('Time [s]');
ylabel('Magnitude');
legend('Car Vertical Velocity (x_1, y)[m/s]', 'Wheel Vertical Velocity (x_2)[m/s]', 'Wheel Vertical Displacement (r)[m]');
title('Simulation of Simplified Suspension System');
grid on;