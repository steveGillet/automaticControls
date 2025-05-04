dataTable = readtable('data.xlsx');
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

figure;
subplot(2,1,1);
semilogx(pResponse(:,3), pResponse(:,1));
grid on;
title('Bode Plot of Empirical Data'); % Title for the entire figure
ylabel('Magnitude (dB)'); % Label for the first subplot
xlim([0.5, 50]);
subplot(2,1,2);
semilogx(pResponse(:,3), pResponse(:,2));
grid on;
xlabel('Frequency (rad/s)'); % Label for the x-axis
ylabel('Phase (deg)'); % Label for the second subplot
xlim([0.5, 50]);

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
figure;

% Magnitude plot
subplot(2,1,1);
semilogx(pResponse(:,3), pResponse(:,1), 'b', 'LineWidth', 1.5); % Empirical data in blue
hold on;
semilogx(w, 20*log10(mag), 'r--', 'LineWidth', 1.5); % Transfer function in red dashed
grid on;
ylabel('Magnitude (dB)');
title('Comparison of Empirical Data and Transfer Function Bode Plot');
legend('Empirical Data', 'Transfer Function');
xlim([0.5, 50]);

% Phase plot
subplot(2,1,2);
semilogx(pResponse(:,3), pResponse(:,2), 'b', 'LineWidth', 1.5); % Empirical data in blue
hold on;
semilogx(w, phase, 'r--', 'LineWidth', 1.5); % Transfer function in red dashed
grid on;
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
legend('Empirical Data', 'Transfer Function');
xlim([0.5, 50]);

[p,z] = pzmap(estTF);
figure;
nyquist(estTF);
xlim([-6, 4]);
ylim([-5, 5]);
disp(p);
disp(z);

figure;
scatter(real(p), imag(p), 'x', 'MarkerEdgeColor', 'r', 'LineWidth', 2); % Poles in red 'x'
hold on;
scatter(real(z), imag(z), 'o', 'MarkerEdgeColor', 'b', 'LineWidth', 2); % Zeros in blue 'o'

xline(0, '--k', 'LineWidth', 1); % Vertical line (imaginary axis)
yline(0, '--k', 'LineWidth', 1); % Horizontal line (real axis)

xlim([-10, 10]); % Adjust as needed to fit your data
ylim([-10, 10]); % Adjust as needed to fit your data

grid on;
xlabel('Real Axis');
ylabel('Imaginary Axis');
title('Manual Pole-Zero Plot');
legend('Poles', 'Zeros');
hold off;

% Convert empirical data to Nyquist plot
realPart = K*10.^(pResponse(:,1) / 20) .* cosd(pResponse(:,2)); % Real part
imagPart = K*10.^(pResponse(:,1) / 20) .* sind(pResponse(:,2)); % Imaginary part

% Plot Nyquist plot of empirical data
figure;
plot(realPart, imagPart, 'b', 'LineWidth', 1.5); % Empirical data in blue
hold on;
plot(realPart, -imagPart, 'b--', 'LineWidth', 1.5); % Mirror for negative frequencies
plot(-1, 0, 'rx', 'MarkerSize', 10, 'LineWidth', 2); % Red "x" at -1
grid on;
xlabel('Real Part');
ylabel('Imaginary Part');
title('Nyquist Diagram of Empirical Data');
legend('Empirical Data', 'Mirror (Negative Frequencies)');
axis equal;

K = 60;
a = 1;
b = 10;
leadNumC = [1 a];
leadDenC = [1 b];
numC = [1 zetaP*omegaP omegaP^2];
denC = [1 zetaZ*omegaZ omegaZ^2];
numC = conv(numC, leadNumC);
denC = conv(denC, leadDenC);
C = tf(K * numC, denC);
closedTF = C*estTF/(1 + C*estTF);
figure;
margin(closedTF);
grid on;
title('Closed Loop Bode Plot with Phase Margin');
figure;
pzmap(closedTF);

% G = frd(10.^(pResponse(:,1)/20) .* exp(1j * deg2rad(pResponse(:,2))), pResponse(:,3));
% numHighRes = [1 0.05*38.6 38.6^2];
% denHighRes = [1 0.1*30.7 30.7^2];
% numC = conv(numC, numHighRes);
% denC = conv(denC, denHighRes);
% C = tf(K * numC, denC);
% empiricalClosedTF = G * C/(1 + G * C);
% figure;
% margin(empiricalClosedTF);

% Define colors
sunset_orange = [1, 0.5765, 0.1608]; % RGB [255, 147, 41]
contrasting_red = [0.7843, 0.1373, 0.1373]; % RGB [200, 35, 35]
torque_gray = [0.3, 0.3, 0.3]; % Dark gray for torque signal

% Plot step response
figure;
subplot(2,1,1);
plot(out.simout1.Time, out.simout1.Data, 'Color', sunset_orange, 'LineWidth', 1.5, 'DisplayName', 'Input r(t)');
hold on;
plot(out.simout.Time, out.simout.Data, 'Color', contrasting_red, 'LineWidth', 1.5, 'DisplayName', 'Output y(t)');
grid on;
xlabel('Time (s)');
ylabel('Output (rad/s)');
title('Step Response (0.5 rad)');
legend('reference r(t)', 'output y(t)');

% Plot torque input signal for step
subplot(2,1,2);
plot(out.simout2.time, out.simout2.Data, 'r', 'LineWidth', 1.5);
hold on;
% Add limiting values (assume ±67 mNm)
yline(67, 'k--', 'Label', 'Upper Limit', 'LineWidth', 1);
yline(-67, 'k--', 'Label', 'Lower Limit', 'LineWidth', 1);
grid on;
xlabel('Time (s)');
ylabel('Torque Input Signal (mNm)');
title('Torque Input Signal for 0.5 Radian Step Input');

% New Model Flipped Resonances
K = 1;
zetaP = 0.05;
zetaZ = 0.1;
omegaP = 4.60118;
omegaZ = 8.34686;

num = K*[1 zetaZ*omegaZ omegaZ^2];
den = [1 zetaP*omegaP omegaP^2 0];

pole = [1 0.65];
den = conv(den, pole);

flippedTF = tf(num,den);
figure;
bode(flippedTF);
xlim([0.5, 50]);
grid on;
title('Bode Plot of Flipped Resonances');

K = 60;
a = 1;
b = 10;
leadNumC = [1 a];
leadDenC = [1 b];
numC = [1 zetaZ*omegaZ omegaZ^2];
denC = [1 zetaP*omegaP omegaP^2];
numC = conv(numC, leadNumC);
denC = conv(denC, leadDenC);
C = tf(K * numC, denC);
closedTF = C*flippedTF/(1 + C*flippedTF);
figure;
margin(closedTF);
grid on;
title('Flipped Closed Loop Bode Plot with Phase Margin');