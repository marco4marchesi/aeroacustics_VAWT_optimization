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
user_settings;



addpath(pwd)
matlab_graphics;



%% simulation data
omega = 21.33; % [
f = omega/(2*pi); %[Hz]
L = 1.5; % [m] diameter
U = 8; % [m/s] , asymptotic flow velocity
St = f*L/U;
Nb = 3;
%% post process - which script do you want to run?
testcase_folder = 'RRF_1p_MACH_0';
flow_RRF_1_nominal = acousticPostProcess([simulationsFolderPath,testcase_folder], 2e-4);

testcase_folder = 'RRF_3p_MACH';
flow_RRF_3_nominal = acousticPostProcess([simulationsFolderPath,testcase_folder], 2e-4);

testcase_folder = 'unsteady_3_profili_nominal';
flow_unst_3_nominal = acousticPostProcess([simulationsFolderPath,testcase_folder], 2e-4);

testcase_folder = 'unsteady_3_profili_deform2';
flow_unst_3_deform2 = acousticPostProcess([simulationsFolderPath,testcase_folder], 1e-4);

testcase_folder = 'unsteady_3_profili_deform6';
flow_unst_3_deform6 = acousticPostProcess([simulationsFolderPath,testcase_folder], 1e-4);

% testcase_folder = 'unsteady_1_profilo_nominal';
% flow_unst_1_nominal = acousticPostProcess([simulationsFolderPath,testcase_folder], 2e-4);

testcase_folder = 'unsteady_1_profilo_deform6';
flow_unst_1_deform6 = acousticPostProcess([simulationsFolderPath,testcase_folder], 2e-4);

%% plots
pixel_x = 900;
pixel_y = 800;


%% 5 metri
% p rms (total)
fig_polar_directivity_p = figure('Position',[100,100,pixel_x,pixel_y]);
polarplot(flow_unst_3_nominal.theta_obs{1},flow_unst_3_nominal.SPL_p{1},'r-*','DisplayName','3p Nominal');
hold on;
polarplot(flow_RRF_1_nominal.theta_obs{1},flow_RRF_1_nominal.SPL_p{1},'g-s','DisplayName','1p RRF nominal');
polarplot(flow_RRF_3_nominal.theta_obs{1},flow_RRF_3_nominal.SPL_p{1},'b-d','DisplayName','3p RRF nominal');

% polarplot(flow_unst_3_deform2.theta_obs{1},flow_unst_3_deform2.SPL_p{1},'g-*','DisplayName','3p Deform 2');
% polarplot(flow_unst_3_deform6.theta_obs{1},flow_unst_3_deform6.SPL_p{1},'b-*','DisplayName','3p Deform 6');
legend
rlim([62,72])
title("Unsteady $SPL$ p rms")
exportgraphics(fig_polar_directivity_p,imagesPath+"11_directivity_unst_vs_rrf.png")

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


%% pressure oscillations vs load oscillations
figure('Position',[100,100,pixel_x,pixel_y])
hold on
% plot(flow_unst_3_nominal.time(2:end),flow_unst_3_nominal.dCL/1000,'--');
% plot(flow_unst_3_nominal.time(2:end),flow_unst_3_nominal.dCD/1000,'--');

plot((flow_unst_3_nominal.pp001r10(:,1)-flow_unst_3_nominal.pp001r10(1,1)),flow_unst_3_nominal.pp001r10(:,5));
plot(flow_unst_3_nominal.time(2:end)-flow_unst_3_nominal.time(2),flow_unst_3_nominal.dCF/1000);
legend('Loading','dCF')
%% wtf is this
figure('Position',[100,100,pixel_x,pixel_y])
hold on
plot(flow_unst_3_nominal.time,flow_unst_3_nominal.pp001r10(:,4));
plot(flow_unst_3_nominal.time,flow_unst_3_nominal.equivalent_speed/100000);
%%
% figure('Position',[100,100,pixel_x,pixel_y])
% for ii=1:120
%     edf(ii) = diff(flow_unst_1_deform6.(strcat('pp',num2str(ii,'%03d'),'r10')){1,2}(1:2,1));
% end
% polarplot([0:3:357]*pi/180,edf)
% rlim([1.98e-4,1.988e-4])

%% FFT analysis
audible_freq = [20,20000];%Hz, minimum and maximum audible frequencies
x1 = [0.01,audible_freq(1)];
y = [0;0.1];
x2 = [audible_freq(2),10000000];
pa1X = [x1, flip(x1)]';
pa2X = [x2, flip(x2)]';
paY = [y(1),y(1),y(2),y(2)]';
alpha_val = 0.3;


fig_fft = figure('Position',[100,100,800,600]);
semilogx(flow_unst_3_nominal.pp001r5_fft_f,flow_unst_3_nominal.pp001r5_fft,'DisplayName','Unsteady')
hold on;
semilogx(flow_RRF_3_nominal.pp031r5_fft_f,flow_RRF_3_nominal.pp031r5_fft,'DisplayName','RRF')
a = patch(pa1X,paY,'r','DisplayName','Inaudible frequencies');
b = patch(pa2X,paY,'r','HandleVisibility','off');
alpha(a,alpha_val);
alpha(b,alpha_val);
legend('Location','northeast')
title('')
xlabel('f [Hz]')
ylabel('p [Pa]')
xline(f,'-.',num2str(f),'DisplayName','1/T','LineWidth',2)
xline(Nb*f,'--',num2str(Nb*f),'DisplayName','Nb/T','LineWidth',2)
xlim([1,10000])
exportgraphics(fig_fft,imagesPath+"fft_acoustic_pressure.png")

% STROUHAL
fig_Strouhal = figure('Position',[100,100,pixel_x,pixel_y]);
semilogx(flow_unst_3_nominal.pp001r5_PSD_St,flow_unst_3_nominal.pp001r5_PSD,'DisplayName','Nominal at 0')
hold on;
semilogx(flow_unst_3_nominal.pp031r5_PSD_St,flow_unst_3_nominal.pp031r5_PSD,'DisplayName','Nominal at 90')
legend
title('')
xlabel('St [-]')
ylabel('p [Pa]')
xline(Nb*St,'--',num2str(Nb*St),'DisplayName','Nb/T','LineWidth',2)
% exportgraphics(fig_Strouhal,imagesPath+"Strouhal_acoustic_pressure.png")
