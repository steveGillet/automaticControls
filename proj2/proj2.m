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

desiredPoles = [-4.44 + 4.44i, -4.44 - 4.44i, -10, -12];
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