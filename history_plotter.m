clear; clc; close all;

folder = "\\wsl.localhost\Ubuntu-20.04\home\marco\aeroProject\";
filename = "history_A.csv";
history = table2array(readtable(folder+filename));

RMSrho_Index = 2;
CauchyCMz_Index = 12;
CMz_Index = 11;
CFL_Index = 7;

figure('Position',[100,100,600,400])
semilogy(history(:,1),history(:,RMSrho_Index))
grid on;
title('RMS \rho')

figure('Position',[100,100,600,400])
semilogy(history(:,1),history(:,CauchyCMz_Index))
grid on;
title('Cauchy CMz')

figure('Position',[100,100,600,400])
plot(history(:,1),history(:,CMz_Index))
grid on;
title('CMz')

figure('Position',[100,100,600,400])
plot(history(:,1),history(:,CFL_Index))
grid on;
title('CFL')
