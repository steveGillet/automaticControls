La = 0.010; % [H]
Ra = 1.0; % [Ohm]
J1 = 0.1; % [kg m^2]
J2 = 1.0; % [kg m^2]
k = 10; % [Nm/rad]
b = 0.01; % [Nms/rad]
Kt = 1.0; % [Nm/A]
Ke = 1.0; % [Vs/rad]
B = 0.01; % [Nms/rad]

A = [-Ra/La 0 -Ke/La 0 0; 0 0 1 0 0; Kt/J1 -k/J1 (-B-b)/J1 k/J1 b/J1; 0 0 0 0 1; 0 k/J2 b/J2 -k/J2 -b/J2]; 
B = [1/La; 0; 0; 0; 0];       
C = [0 1 0 0 0];        
D = [0];  

rotorSS = ss(A, B, C, D);

disp(eig(rotorSS));

t = 0:0.01:10; 
Va = zeros(size(t)); 
Va(t >= 1) = 1;
u = [Va];

[y, t, x] = lsim(rotorSS, u, t);

% Plot the input and output
figure;
plot(t, u, 'g-', t, y, 'r-');
xlabel('Time [s]');
ylabel('Magnitude');
legend('Voltage Armature (u) [V]', 'Shaft Angle (Theta1) [rad]');
title('Simulation of Electric Motor Driving Flexible Load');
grid on;