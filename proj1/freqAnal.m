dataTable = readtable('data.xlsx');
data = table2array(dataTable);

n = size(data,1);
pResponse = zeros(n,3);

for i = 1:n
    pResponse(i,1) = 20*log10(data(i,2)) - 20*log10(data(i,1));
    pResponse(i,2) = data(i,3)*180/pi - 90;
    pResponse(i,3) = data(i,1)*2*pi;
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