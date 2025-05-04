dataTable = readtable('data.xlsx', 'VariableNamingRule', 'preserve');
data = table2array(dataTable);

[uniqueFreq, ia] = unique(data(:,1));
n = size(uniqueFreq,1);
pResponse = zeros(n,3);

for i = 1:n
    idx = ia(i);
    pResponse(i,1) = 20*log10(data(idx,2)) - 20*log10(data(i,1));
    pResponse(i,2) = data(idx,3)*180/pi - 90;
    pResponse(i,3) = data(idx,1)*2*pi;
end

K = 1;
% K = -0.05;
zetaZ = 0.05;
zetaP = 0.1;
omegaZ = 4.60118;
omegaP = 8.34686;

num = K*[1 zetaZ*omegaZ omegaZ^2];
den = [1 zetaP*omegaP omegaP^2 0];
pole = [1 0.65];
den = conv(den, pole);
% pole = [1 1];
% den = conv(den, pole);
estTF = tf(num,den);

% Generate Bode plot data for the transfer function
[mag, phase, w] = bode(estTF);
mag = squeeze(mag); % Convert to 1D array
phase = squeeze(phase); % Convert to 1D array

% Plot empirical data and transfer function Bode plot on the same figure
clf; figure(1); 

% Magnitude plot
subplot(2,1,1);
semilogx(pResponse(:,3), pResponse(:,1), 'b', 'LineWidth', 1.5); % Empirical data in blue
hold on;
semilogx(w, 20*log10(mag), 'r--', 'LineWidth', 1.5); % Transfer function in red dashed
grid on;
ylabel('Magnitude (dB)');
title('Comparison of Empirical Data and Analytical Transfer Function Bode Plot');
legend('Empirical Data', 'Analytical Transfer Function');
xlim([0.5, 50]);

% Phase plot
subplot(2,1,2);
semilogx(pResponse(:,3), pResponse(:,2), 'b', 'LineWidth', 1.5); % Empirical data in blue
hold on;
semilogx(w, phase, 'r--', 'LineWidth', 1.5); % Transfer function in red dashed
grid on;
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
legend('Empirical Data', 'Analytical Transfer Function');
xlim([0.5, 50]);

minSS = ss(estTF);
minSS = minreal(minSS);
disp(minSS.C);

[A,B,C,D] = tf2ss(num, den);
disp(A);
disp(B);
disp(C);
disp(D);

disp(eig(A));

desiredPoles = [-0.615 + 4.5997i, -0.615 - 4.5997i, -10, -12];
K = place(A, B, desiredPoles);
F = inv(C*inv(-A+B*K)*B);

sysCL = ss(A-B*K, B*F, C, D);
figure(2);
bode(sysCL);
title('Closed Loop Tracking Frequency Response');
grid on;

sysControl = ss(A-B*K, B*F, -K, 0);
figure(3);
bode(sysControl);
title('Control Effort Bode Plot');

sysLG = ss(A,B,K,0);
figure(4);
bode(sysLG);
title('Loop Gain Bode Plot');
grid on;

figure(5);
nyquist(sysLG);
title('Loop Gain Nyquist Plot');

figure(6);
margin(sysLG);

[p,z] = pzmap(estTF);
disp(p);
disp(z);

figure(7);
margin(sysControl);
disp(K);

% Define sunset colors
sunsetOrange = [1, 0.4353, 0.3804]; % For reference r
sunsetPurple = [0.4196, 0.3569, 0.5843]; % For output y
sunsetPink = [1, 0.4118, 0.7059]; % For control input u

% Extract time and signals
time = out.r.Time; % Time vector (seconds)
r = out.r.Data; % Reference input r (rad)
y = out.y.Data; % Plant output y (rad)
u = out.u.Data; % Control input u (mN-m)

% Plot 1: Reference r and Plant Output y vs. Time on the same plot
figure(8);
subplot(2,1,1);
plot(time, r, 'Color', sunsetOrange, 'LineWidth', 1.5, 'DisplayName', 'Reference r(t)');
hold on;
plot(time, y, 'Color', sunsetPurple, 'LineWidth', 1.5, 'DisplayName', 'Output y(t)');
hold off;
grid on;
title('Reference and Plant Output Response');
xlabel('Time (s)');
ylabel('Angular Deflection (rad)');
legend('show');

% Plot 2: Control Input u vs. Time
subplot(2,1,2);
plot(time, u, 'Color', sunsetPink, 'LineWidth', 1.5);
grid on;
title('Control Input u(t)');
xlabel('Time (s)');
ylabel('Torque (mN-m)');

% Adjust plot layout
sgtitle('Simulation Results for Sine Reference Input [0.1 rad at 0.1 Hz]');

scalingFactor = 5;
desiredObserverPoles = [scalingFactor*-0.615 + 4.5997i, scalingFactor*-0.615 - 4.5997i, scalingFactor*-10, scalingFactor*-12];

disp(desiredObserverPoles);

L = place(A', C', desiredObserverPoles)';

obsSys = ss(A-B*K-L*C, L, -K, 0);
analG = ss(A, B, C, D);

obsLGsys = -obsSys*analG;
figure(9);
margin(obsLGsys);
title('Analytical Observer Loop Gain Bode Plot');
grid on;
xlim([0.5, 50]);

figure(10);
nyquist(obsLGsys);
title('Analytical Observer Loop Gain Nyquist Plot');

empG = frd(10.^(pResponse(:,1)/20) .* exp(1j * deg2rad(pResponse(:,2))), pResponse(:,3));

empObsLGsys = -obsSys*empG;
figure(11);
margin(empObsLGsys);
title('Empirical Observer Loop Gain Bode Plot');
grid on;
xlim([0.5, 50]);

figure(12);
nyquist(empObsLGsys);
title('Empirical Observer Loop Gain Nyquist Plot');
