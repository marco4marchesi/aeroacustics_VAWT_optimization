%{

test script for aeroacoustic fwh post process

structure of the directories:
- results saved in the "results" directory
- surface_flow.dat in the "surfaceFiles" directory
- matlab scripts in the "postProcess" folder: in this folder you will have
all your matlab scripts

%}

%% init - do not modify
clear;  clc; close all;
main_folder = pwd;
user = 'marco';
imagesPath = "C:\Users\marco\OneDrive - Politecnico di Milano\MAGISTRALE\QuartoSemestre\Aeroacoustics\aeroacoustic project\varie";
addpath(pwd)
matlab_graphics;


%% simulation data
omega = 21.33;
f = omega/(2*pi);

%% post process - which script do you want to run?
testcase_folder = 'unsteady_3_profili_deform6';
su2_folder = ['E:\UNI - fisso\aeroacustica\',testcase_folder];
flow_unst_3_deform6 = acousticPostProcess(su2_folder, 2e-4, 1472);

% testcase_folder = 'unsteady_3_profili_deform6';
% su2_folder = ['\\wsl.localhost\ubuntu-20.04\home\',user,'\',testcase_folder];
% flow_unst_3_deform6 = acousticPostProcess(su2_folder, 2e-4, 1472);
% 
% testcase_folder = 'unsteady_3_profili_nominal';
% su2_folder = ['\\wsl.localhost\ubuntu-20.04\home\',user,'\',testcase_folder];
% flow_unst_3_nominal = acousticPostProcess(su2_folder, 2e-4, 1472);


%% plots
pixel_x = 600;
pixel_y = 500;


%% 5 metri
fig_polar_directivity_p = figure('Position',[100,100,pixel_x,pixel_y]);
% polarplot(flow_unst_3_nominal.theta_obs{1},flow_unst_3_nominal.SPL_p{1},'r-*','DisplayName','Nominal');
% hold on;
% polarplot(flow_unst_3_deform1.theta_obs{1},flow_unst_3_deform1.SPL_p{1},'g-*','DisplayName','Deform 1');
polarplot(flow_unst_3_deform6.theta_obs{1},flow_unst_3_deform6.SPL_p{1},'b-*','DisplayName','Deform 6');
legend
rlim([50,75])
title("Unsteady $SPL$ p rms")
% exportgraphics(fig_polar_directivity_p,imagesPath+"/directivity_p.png")
pointMatrix.X = [flow_unst_3_deform6.points{1}(:,1),flow_unst_3_deform6.points{2}(:,1),flow_unst_3_deform6.points{3}(:,1)];
pointMatrix.Y = [flow_unst_3_deform6.points{1}(:,2),flow_unst_3_deform6.points{2}(:,2),flow_unst_3_deform6.points{3}(:,2)];
pointMatrix.Z = [flow_unst_3_deform6.SPL_p{1}',flow_unst_3_deform6.SPL_p{2}',flow_unst_3_deform6.SPL_p{3}'];

fig_3D_directivity = figure('Position',[100,100,pixel_x,pixel_y]);
surf(pointMatrix.X,pointMatrix.Y, pointMatrix.Z)
axis equal
xlabel('x [m]')
ylabel('y [m]')
zlabel('SPL(p) [dB]')
colorbar
% fig_polar_directivity_thickness = figure('Position',[100,100,pixel_x,pixel_y]);
% polarplot(flow_unst_3_nominal.theta_obs{1},flow_unst_3_nominal.SPL_thickness{1},'r-*','DisplayName','Nominal');
% hold on;
% polarplot(flow_unst_3_deform1.theta_obs{1},flow_unst_3_deform1.SPL_thickness{1},'g-*','DisplayName','Deform 1');
% polarplot(flow_unst_3_deform6.theta_obs{1},flow_unst_3_deform6.SPL_thickness{1},'b-*','DisplayName','Deform 6');
% thetalim([0 360])
% rlim([20 32])
% legend
% title("$SPL$ thickness rms")
% % % exportgraphics(fig_polar_directivity_thickness,imagesPath+"/directivity_thickness.png")
% 
% 
% 
% fig_polar_directivity_loading = figure('Position',[100,100,pixel_x,pixel_y]);
% polarplot(flow_unst_3_nominal.theta_obs{1},flow_unst_3_nominal.SPL_loading{1},'r-*','DisplayName','Nominal');
% hold on;
% polarplot(flow_unst_3_deform1.theta_obs{2},flow_unst_3_deform1.SPL_loading{1},'g-*','DisplayName','Deform 1');
% polarplot(flow_unst_3_deform6.theta_obs{3},flow_unst_3_deform6.SPL_loading{1},'b-*','DisplayName','Deform 6');
% thetalim([0 360])
% rlim([50 80])
% legend
% title("$SPL$ loading rms")
% exportgraphics(fig_polar_directivity_loading,imagesPath+"/directivity_loading.png")







% fig_acoustic_pressure = figure('Position',[100,100,pixel_x,pixel_y]);
% plot(flow_unst_3_deform1.pp001(:,1),flow_unst_3_deform1.pp001(:,3),'DisplayName','unst def1')
% hold on;
% plot(flow_unst_3_deform6.pp001(:,1),flow_unst_3_deform6.pp001(:,3),'DisplayName','unst def6')
% 
% % plot(flow_rrf_single.pp001(:,1),flow_rrf_single.pp001(:,3),'MarkerSize',8,'DisplayName','RRF 1')
% % plot(flow_unsteady.pp001(:,1),flow_unsteady.pp001(:,3),'MarkerSize',8,'DisplayName','Unsteady')
% legend
% xlabel('Time [s]')
% ylabel('Acoustic pressure [Pa]')
% exportgraphics(fig_acoustic_pressure,imagesPath+"/acoustic_pressure.png")

%% FFT analysis

% [flow_unst_3_deform6.pp001_fft_f,flow_unst_3_deform6.pp001_fft] = fourierSingleSided(1/(2e-4), flow_unst_3_deform6.pp001(last_lap_indexes,3)-mean(flow_unst_3_deform6.pp001(last_lap_indexes,3)));
% fig_fft = figure('Position',[100,100,pixel_x,pixel_y]);
% % plot(flow_unsteady.pp001_fft_f,flow_unsteady.pp001_fft,'DisplayName','Unsteady')
% plot(flow_unst_3_deform1.pp001_fft_f,flow_unst_3_deform1.pp001_fft,'DisplayName','unst def')
% hold on;
% legend
% title('')
% xlabel('f [Hz]')
% ylabel('p [Pa]')
% xline(3*f,'--',num2str(3*f),'DisplayName','3/T','LineWidth',2)
% xline(4*f,'b-.',num2str(4*f),'DisplayName','4/T','LineWidth',2)
% xlim([0,100])
% exportgraphics(fig_fft,imagesPath+"/fft_acoustic_pressure.png")

