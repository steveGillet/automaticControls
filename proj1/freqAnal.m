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
