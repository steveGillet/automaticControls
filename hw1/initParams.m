coeff3=0.000028291; % [Vs^3]
coeff2=0.285890679; % [Vs^2]
coeff1=1.332; % [Vs/rad]
num = [-1];
den = [coeff3, coeff2, coeff1, 0];
format longG;
poles = roots(den);
disp(poles);

figure;
hold on;
plot(real(poles), imag(poles), 'bx', 'MarkerSize', 10, 'LineWidth', 2);
title('Root Locus G 0 to -1');
xlabel('Re');
ylabel('Im');
grid on;

Gp = -10;
Gd = -0.1;
numCL = [-Gd -Gp];
denCL = [coeff3, coeff2, coeff1 - Gd, -Gp];
polesCL = roots(denCL);
disp(polesCL);

plot(real(polesCL), imag(polesCL), 'rx', 'MarkerSize', 10, 'LineWidth', 2);


gScalar = linspace(0, -1, 100);
for i = 1:length(gScalar)
    g = gScalar(i);
    denCL = [coeff3, coeff2, coeff1 - Gd*g, -Gp*g];
    polesCL = roots(denCL);
    plot(real(polesCL), imag(polesCL), 'k.', 'MarkerSize', 1);
    disp(polesCL);
end

legend('Open Loop Poles', 'Closed Loop Poles', 'Root Locus');
hold off;

plot(out.closedIn.time, out.closedIn.Data, 'b-', ...
     out.closedOut.time, out.closedOut.Data, 'r--');

xlabel('Time (s)');
ylabel('Arm Angle (rads)');
legend('Input Theta_R', 'Output Theta_L');
title('-1 g Closed Loop Arm Angle Tracking');
grid on;