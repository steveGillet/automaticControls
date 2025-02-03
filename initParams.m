coeff3=0.000266135; % [Vs^3]
coeff2=2.6893656; % [Vs^2]
coeff1=0.132; % [Vs/rad]
num = [-1];
den = [coeff3, coeff2, coeff1, 0];
format longG;
disp(roots(den));

plot(out.inData.time, out.inData.Data, 'b-', ...
     out.outData.time, out.outData.Data, 'r--');

xlabel('Time (s)');
ylabel('Voltage (V)');
legend('Input Vp', 'Sensor Voltage Vs');
title('Input vs. Output');
grid on;