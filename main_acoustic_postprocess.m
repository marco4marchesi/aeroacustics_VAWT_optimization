%{

test script for aeroacoustic fwh post process

structure of the directories:
- results saved in the "results" directory
- surface_flow.dat in the "surfaceFiles" directory
- matlab scripts in the "postProcess" folder: in this folder you will have
all your matlab scripts

% simulationsFolderPath:
% this changes depending on your pc, I have all the simulations as
% subfolder (test case) in a single mother folder that I put here:
% remember to put the "\" at the end of the string...


%}

%% init - do not modify
clear;  clc; close all;
main_folder = pwd;

% who uses this script? select user
user = 'marco';


switch user
    case 'marco'
        imagesPath = "C:\Users\marco\OneDrive - Politecnico di Milano\MAGISTRALE\QuartoSemestre\Aeroacoustics\aeroacoustic project\varie";
        simulationsFolderPath = 'E:\UNI - fisso\aeroacustica\';
    case 'fra'
        imagesPath = [];
        simulationsFolderPath = [];
    case 'adri'
        imagesPath = [];
        simulationsFolderPath = [];
end

addpath(pwd)
matlab_graphics;



%% simulation data
omega = 21.33;
f = omega/(2*pi);

%% post process - which script do you want to run?
testcase_folder = 'unsteady_3_profili_nominal';
flow_unst_3_nominal = acousticPostProcess([simulationsFolderPath,testcase_folder], 2e-4, 1472);

testcase_folder = 'unsteady_3_profili_deform2';
flow_unst_3_deform2 = acousticPostProcess([simulationsFolderPath,testcase_folder], 1e-4, 2944);

testcase_folder = 'unsteady_3_profili_deform6';
flow_unst_3_deform6 = acousticPostProcess([simulationsFolderPath,testcase_folder], 1e-4, 2944);

%% plots
pixel_x = 600;
pixel_y = 500;


%% 5 metri
% p rms (total)
fig_polar_directivity_p = figure('Position',[100,100,pixel_x,pixel_y]);
polarplot(flow_unst_3_nominal.theta_obs{1},flow_unst_3_nominal.SPL_p{1},'r-*','DisplayName','3p Nominal');
hold on;
polarplot(flow_unst_3_deform2.theta_obs{1},flow_unst_3_deform2.SPL_p{1},'g-*','DisplayName','3p Deform 2');
polarplot(flow_unst_3_deform6.theta_obs{1},flow_unst_3_deform6.SPL_p{1},'b-*','DisplayName','3p Deform 6');
legend
title("Unsteady $SPL$ p rms")
% exportgraphics(fig_polar_directivity_p,imagesPath+"/directivity_p.png")

% thickness
fig_polar_directivity_thickness = figure('Position',[100,100,pixel_x,pixel_y]);
polarplot(flow_unst_3_nominal.theta_obs{1},flow_unst_3_nominal.SPL_thickness{1},'r-*','DisplayName','Nominal');
hold on;
polarplot(flow_unst_3_deform2.theta_obs{1},flow_unst_3_deform2.SPL_thickness{1},'g-*','DisplayName','Deform 2');
polarplot(flow_unst_3_deform6.theta_obs{1},flow_unst_3_deform6.SPL_thickness{1},'b-*','DisplayName','Deform 6');
legend
title("$SPL$ thickness rms")
% exportgraphics(fig_polar_directivity_thickness,imagesPath+"/directivity_thickness.png")
 
% loading
fig_polar_directivity_loading = figure('Position',[100,100,pixel_x,pixel_y]);
polarplot(flow_unst_3_nominal.theta_obs{1},flow_unst_3_nominal.SPL_loading{1},'r-*','DisplayName','Nominal');
hold on;
polarplot(flow_unst_3_deform2.theta_obs{1},flow_unst_3_deform2.SPL_loading{1},'g-*','DisplayName','Deform 2');
polarplot(flow_unst_3_deform6.theta_obs{1},flow_unst_3_deform6.SPL_loading{1},'b-*','DisplayName','Deform 6');
legend
title("$SPL$ loading rms")
% exportgraphics(fig_polar_directivity_loading,imagesPath+"/directivity_loading.png")

% 3D directivity
fig_3D_directivity = figure('Position',[100,100,pixel_x,pixel_y]);
surf(flow_unst_3_deform6.pointMatX,flow_unst_3_deform6.pointMatY, flow_unst_3_deform6.p_rmsMat)
axis equal
xlabel('x [m]')
ylabel('y [m]')
zlabel('SPL(p) [dB]')
colorbar
% exportgraphics(fig_3D_directivity,imagesPath+"/3D_directivity_p.png")



%%
% fig_acoustic_pressure = figure('Position',[100,100,pixel_x,pixel_y]);
% plot(flow_unst_3_deform1.pp001(:,1),flow_unst_3_deform1.pp001(:,3),'DisplayName','unst def1')
% hold on;
% plot(flow_unst_3_deform6.pp001(:,1),flow_unst_3_deform6.pp001(:,3),'DisplayName','unst def6')

% plot(flow_rrf_single.pp001(:,1),flow_rrf_single.pp001(:,3),'MarkerSize',8,'DisplayName','RRF 1')
% plot(flow_unsteady.pp001(:,1),flow_unsteady.pp001(:,3),'MarkerSize',8,'DisplayName','Unsteady')
% exportgraphics(fig_acoustic_pressure,imagesPath+"/acoustic_pressure.png")

%% FFT analysis
% 
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
%%
figure
hold on
% plot(flow_unst_3_nominal.time(2:end),flow_unst_3_nominal.dCL/1000,'--');
% plot(flow_unst_3_nominal.time(2:end),flow_unst_3_nominal.dCD/1000,'--');
plot(flow_unst_3_nominal.time,flow_unst_3_nominal.pp010r10{1,2}(:,5));
plot(flow_unst_3_nominal.time(2:end),flow_unst_3_nominal.dCF/1000);
legend('Acoustic Pressure','dCF')
%%
figure
hold on
plot(flow_unst_3_nominal.time,flow_unst_3_nominal.pp010r10{1,2}(:,4));
plot(flow_unst_3_nominal.time,flow_unst_3_nominal.equivalent_speed/100000);
%%
for ii=1:120
    edf(ii) = diff(flow_unst_3_nominal.(strcat('pp',num2str(ii,'%03d'),'r10')){1,2}(1:2,1));
end
polarplot([0:3:357]*pi/180,edf)
rlim([3.94e-4,3.9451e-4])