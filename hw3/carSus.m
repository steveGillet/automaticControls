b = 100000;
m1 = 10; % [kg]
m2 = 350; % [kg]
kw = 500000; % [Nm]
ks = 10000; % [Nm]


A = [-b/m2 -1/m2 b/m2 0; ks 0 -ks 0; 0 0 0 -1/m1; 0 0 kw 0]; 
B = [0; 0; 0; -kw];       
C = [1 0 0 0];        
D = [0];  

carSusSS = ss(A, B, C, D);
carSusTF = tf(carSusSS);
disp(carSusTF);

plot(wheelDisplacement.time, wheelDisplacement.Data, 'b-', ...
     carDisplacement.time, carDisplacement.Data, 'r-');

xlabel('Time [s]');
ylabel('Displacement [m]');
legend('Wheel Displacement', 'Car Displacement');
title(sprintf('b = %d, Wheel vs Car Displacement',b));
grid on;