dataTable = readtable('data.xlsx');
data = table2array(dataTable);

disp(data(1,1));
n = size(data,1);
pResponse = zeros(n,2);

for i = 1:n
    pResponse(i,1) = 20*log10(data(i,2)) - 20*log10(data(i,1));
    pResponse(i,2) = data(i,3)*180/pi - 90;
end

figure;
subplot(2,1,1);
semilogx(data(:,1), pResponse(:,1));
subplot(2,1,2);
semilogx(data(:,1), pResponse(:,2));

K = 0.20;
zetaZ = 0.025;
zetaP = 0.05;
omegaZ = 0.7323;
omegaP = 1.32844;

num = K*[1 zetaZ*omegaZ omegaZ^2];
den = [1 zetaP*omegaP omegaP^2 0];
pole = [2 1];
den = conv(den, pole);
pole = [2 1];
den = conv(den, pole);
estTF = tf(num,den);
figure;
bode(estTF);
disp(estTF);