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

%% post process - EXTRACT DATA
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

testcase_folder = 'unsteady_1_profilo_nominal';
flow_unst_1_nominal = acousticPostProcess([simulationsFolderPath,testcase_folder], 4e-4);

testcase_folder = 'unsteady_1_profilo_deform6';
flow_unst_1_deform6 = acousticPostProcess([simulationsFolderPath,testcase_folder], 2e-4);

testcase_folder = 'unsteady_nominal_hybrid';
flow_unst_hybrid = acousticPostProcess([simulationsFolderPath,testcase_folder], 2e-4);

%% plots
factor = 80; % define the dimensions
pixel_x = 11*factor; % define the proportions
pixel_y = 8*factor;


%% rrf vs unsteady
% 1p rrf vs unsteady p rms (total)
fig_polar_directivity_p_1p = figure('Position',[100,100,pixel_x,pixel_y]);
polarplot(flow_unst_1_nominal.theta_obs{1},flow_unst_1_nominal.SPL_p{1},'r-*','DisplayName','1p Nominal');
hold on;
polarplot(flow_RRF_1_nominal.theta_obs{1},flow_RRF_1_nominal.SPL_p{1},'g-s','DisplayName','1p RRF nominal');
legend('Location','northeast')
rlim([60,73])
title("Unsteady $SPL$ 1p p rms")
% exportgraphics(fig_polar_directivity_p_1p,imagesPath+"11_directivity_unst_vs_rrf_1p.emf")

% 3p rrf vs unsteady prms
fig_polar_directivity_p_3p = figure('Position',[100,100,pixel_x,pixel_y]);
polarplot(flow_unst_3_nominal.theta_obs{1},flow_unst_3_nominal.SPL_p{1},'r-*','DisplayName','3p Nominal');
hold on;
polarplot(flow_RRF_3_nominal.theta_obs{1},flow_RRF_3_nominal.SPL_p{1},'g-s','DisplayName','3p RRF nominal');
legend('Location','northeast')
rlim([60,73])
title("Unsteady $SPL$ 3p p rms")
% exportgraphics(fig_polar_directivity_p_3p,imagesPath+"11_directivity_unst_vs_rrf_3p.emf")

%% DEFORM 2
% p rms
fig_polar_directivity_p_deform2 = figure('Position',[100,100,pixel_x,pixel_y]);
polarplot(flow_unst_3_nominal.theta_obs{1},flow_unst_3_nominal.SPL_p{1},'r-*','DisplayName','Nominal');
hold on;
polarplot(flow_unst_3_deform2.theta_obs{1},flow_unst_3_deform2.SPL_p{1},'g-s','DisplayName','Deform 2');
legend('Location','northeast')
rlim([60,73])
title("Unsteady $SPL$ 3p p rms")
% exportgraphics(fig_polar_directivity_p_deform2,imagesPath+"21_directivity_nom_vs_deform2.emf")

% thickness
fig_polar_directivity_thickness_deform2 = figure('Position',[100,100,pixel_x,pixel_y]);
polarplot(flow_unst_3_nominal.theta_obs{1},flow_unst_3_nominal.SPL_thickness{1},'r-*','DisplayName','Nominal');
hold on;
polarplot(flow_unst_3_deform2.theta_obs{1},flow_unst_3_deform2.SPL_thickness{1},'g-*','DisplayName','Deform 2');
legend('Location','northeast')
title("$SPL$ thickness rms")
rlim([10,30])
% exportgraphics(fig_polar_directivity_thickness_deform2,imagesPath+"/21_directivity_thickness_nom_vs_deform2.emf")
 
% loading
fig_polar_directivity_loading_deform2 = figure('Position',[100,100,pixel_x,pixel_y]);
polarplot(flow_unst_3_nominal.theta_obs{1},flow_unst_3_nominal.SPL_loading{1},'r-*','DisplayName','Nominal');
hold on;
polarplot(flow_unst_3_deform2.theta_obs{1},flow_unst_3_deform2.SPL_loading{1},'g-*','DisplayName','Deform 2');
legend('Location','northeast')
title("$SPL$ loading rms")
rlim([60,80])
% exportgraphics(fig_polar_directivity_loading_deform2,imagesPath+"/21_directivity_loading_nom_vs_deform2.emf")

%% DEFORM 6
% p rms
fig_polar_directivity_p_deform6 = figure('Position',[100,100,pixel_x,pixel_y]);
polarplot(flow_unst_3_nominal.theta_obs{1},flow_unst_3_nominal.SPL_p{1},'r-*','DisplayName','Nominal');
hold on;
polarplot(flow_unst_3_deform6.theta_obs{1},flow_unst_3_deform6.SPL_p{1},'g-s','DisplayName','Deform 6');
legend('Location','northeast')
rlim([60,88])
% title("Unsteady $SPL$ 3p p rms")
% exportgraphics(fig_polar_directivity_p_deform6,imagesPath+"23_directivity_nom_vs_deform6.emf")

% thickness deform 6
fig_polar_directivity_thickness_deform6 = figure('Position',[100,100,pixel_x,pixel_y]);
polarplot(flow_unst_3_nominal.theta_obs{1},flow_unst_3_nominal.SPL_thickness{1},'r-*','DisplayName','Nominal');
hold on;
polarplot(flow_unst_3_deform6.theta_obs{1},flow_unst_3_deform6.SPL_thickness{1},'g-*','DisplayName','Deform 6');
legend('Location','northeast')
title("$SPL$ thickness rms")
rlim([10,50])
% exportgraphics(fig_polar_directivity_thickness_deform6,imagesPath+"/23_directivity_thickness_nom_vs_deform6.emf")
 
% loading deform 6
fig_polar_directivity_loading_deform6 = figure('Position',[100,100,pixel_x,pixel_y]);
polarplot(flow_unst_3_nominal.theta_obs{1},flow_unst_3_nominal.SPL_loading{1},'r-*','DisplayName','Nominal');
hold on;
polarplot(flow_unst_3_deform6.theta_obs{1},flow_unst_3_deform6.SPL_loading{1},'g-*','DisplayName','Deform 6');
legend('Location','northeast')
title("$SPL$ loading rms")
rlim([50,90])
% exportgraphics(fig_polar_directivity_loading_deform6,imagesPath+"/23_directivity_loading_nom_vs_deform6.emf")


%% 1 profile vs 3 profile computations only on 1
% p rms
fig_polar_directivity_hyb_vs_1p = figure('Position',[100,100,pixel_x,pixel_y]);
polarplot(flow_unst_1_nominal.theta_obs{1},flow_unst_1_nominal.SPL_p{1},'r-*','DisplayName','1p');
hold on;
polarplot(flow_unst_hybrid.theta_obs{1},flow_unst_hybrid.SPL_p{1},'g-s','DisplayName','3p 1surf');
legend('Location','northeast')
rlim([60,88])
% title("Unsteady $SPL$ 3p p rms")
% exportgraphics(fig_polar_directivity_hyb_vs_1,imagesPath+"XX_directivity_hyb_vs_1p.emf")

% thickness deform 6
fig_polar_directivity_thickness_hyb_vs_1p= figure('Position',[100,100,pixel_x,pixel_y]);
polarplot(flow_unst_1_nominal.theta_obs{1},flow_unst_1_nominal.SPL_thickness{1},'r-*','DisplayName','1p');
hold on;
polarplot(flow_unst_hybrid.theta_obs{1},flow_unst_hybrid.SPL_thickness{1},'g-*','DisplayName','3p 1surf');
legend('Location','northeast')
title("$SPL$ thickness rms")
rlim([10,50])
% exportgraphics(fig_polar_directivity_thickness_hyb_vs_1p,imagesPath+"/XX_directivity_thickness_hyb_vs_1p.emf")
 
% loading deform 6
fig_polar_directivity_loading_hyb_vs_1p = figure('Position',[100,100,pixel_x,pixel_y]);
polarplot(flow_unst_1_nominal.theta_obs{1},flow_unst_1_nominal.SPL_loading{1},'r-*','DisplayName','1p');
hold on;
polarplot(flow_unst_hybrid.theta_obs{1},flow_unst_hybrid.SPL_loading{1},'g-*','DisplayName','3p 1surf');
legend('Location','northeast')
title("$SPL$ loading rms")
rlim([50,90])
% exportgraphics(fig_polar_directivity_loading_hyb_vs_1p,imagesPath+"/XX_directivity_loading_nom_vs_hyb_vs_1p.emf")


%% 3D directivity
fig_3D_directivity = figure('Position',[100,100,pixel_x,pixel_y]);
surf(flow_unst_3_deform6.pointMatX,flow_unst_3_deform6.pointMatY, flow_unst_3_deform6.p_rmsMat)
axis equal
xlabel('x [m]')
ylabel('y [m]')
zlabel('SPL(p) [dB]')
colorbar
% exportgraphics(fig_3D_directivity,imagesPath+"/3D_directivity_p.emf")


%% COMPARISON OBSERVER DISTANCES
% p rms
fig_polar_directivity_p_distance = figure('Position',[100,100,pixel_x,pixel_y]);
polarplot(flow_unst_3_nominal.theta_obs{1},flow_unst_3_nominal.SPL_p{1},'r-o','DisplayName','R = 5');
hold on;
polarplot(flow_unst_3_nominal.theta_obs{2},flow_unst_3_nominal.SPL_p{2},'g-s','DisplayName','R = 10');
polarplot(flow_unst_3_nominal.theta_obs{3},flow_unst_3_nominal.SPL_p{3},'b-*','DisplayName','R = 20');
legend('Location','northeast')
rlim([40,80])
% title("Unsteady $SPL$ 3p p rms")
% exportgraphics(fig_polar_directivity_p_distance,imagesPath+"20_directivity_nom_distance.emf")

% loading
fig_polar_directivity_loading_distance = figure('Position',[100,100,pixel_x,pixel_y]);
polarplot(flow_unst_3_nominal.theta_obs{1},flow_unst_3_nominal.SPL_loading{1},'r-o','DisplayName','R = 5');
hold on;
polarplot(flow_unst_3_nominal.theta_obs{2},flow_unst_3_nominal.SPL_loading{2},'g-s','DisplayName','R = 10');
polarplot(flow_unst_3_nominal.theta_obs{3},flow_unst_3_nominal.SPL_loading{3},'b-*','DisplayName','R = 20');
legend('Location','northeast')
rlim([40,80])
% title("Unsteady $SPL$ 3p p rms")
% exportgraphics(fig_polar_directivity_loading_distance,imagesPath+"20_directivity_nom_loading_distance.emf")

% thickness
fig_polar_directivity_thickness_distance = figure('Position',[100,100,pixel_x,pixel_y]);
polarplot(flow_unst_3_nominal.theta_obs{1},flow_unst_3_nominal.SPL_thickness{1},'r-o','DisplayName','R = 5');
hold on;
polarplot(flow_unst_3_nominal.theta_obs{2},flow_unst_3_nominal.SPL_thickness{2},'g-s','DisplayName','R = 10');
polarplot(flow_unst_3_nominal.theta_obs{3},flow_unst_3_nominal.SPL_thickness{3},'b-*','DisplayName','R = 20');
legend('Location','northeast')
rlim([-40,30])
% title("Unsteady $SPL$ 3p p rms")
% exportgraphics(fig_polar_directivity_thickness_distance,imagesPath+"20_directivity_nom_thickness_distance.emf")

%%
% fig_acoustic_pressure = figure('Position',[100,100,pixel_x,pixel_y]);
% plot(flow_unst_3_deform1.pp001(:,1),flow_unst_3_deform1.pp001(:,3),'DisplayName','unst def1')
% hold on;
% plot(flow_unst_3_deform6.pp001(:,1),flow_unst_3_deform6.pp001(:,3),'DisplayName','unst def6')

% plot(flow_rrf_single.pp001(:,1),flow_rrf_single.pp001(:,3),'MarkerSize',8,'DisplayName','RRF 1')
% plot(flow_unsteady.pp001(:,1),flow_unsteady.pp001(:,3),'MarkerSize',8,'DisplayName','Unsteady')
% exportgraphics(fig_acoustic_pressure,imagesPath+"/acoustic_pressure.emf")


%% pressure oscillations vs load oscillations
fig_loading_vs_dLoad = figure('Position',[100,100,pixel_x,pixel_y]);
yyaxis left
plot((flow_unst_3_nominal.pp001r10(:,1)-flow_unst_3_nominal.pp001r10(1,1)),flow_unst_3_nominal.pp001r10(:,5)-mean(flow_unst_3_nominal.pp001r10(:,5)),'DisplayName','Loading at 0°');
ylabel('P'' [Pa]')
yyaxis right
plot(flow_unst_3_nominal.time(2:end)-flow_unst_3_nominal.time(2),(flow_unst_3_nominal.dCF-mean(flow_unst_3_nominal.dCF)),'DisplayName','dCF');
ylabel('$\frac{d}{dt}(C_F)$ [1/s]')
xlabel('Time [s]')
legend
% exportgraphics(fig_loading_vs_dLoad,imagesPath+"19_fig_loading_vs_dLoad.emf")


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
y = [1e-6;0.3];
x2 = [audible_freq(2),10000000];
pa1X = [x1, flip(x1)]';
pa2X = [x2, flip(x2)]';
paY = [y(1),y(1),y(2),y(2)]';
alpha_val = 0.1;



fig_fft = figure('Position',[100,100,pixel_x,pixel_y]);
% subplot(2,1,1)
loglog(flow_unst_3_nominal.pp001r5_fft_f/omega*2*pi,flow_unst_3_nominal.pp001r5_fft,'DisplayName','Unsteady $p_{rms}$')
hold on;
loglog(flow_RRF_3_nominal.pp001r5_fft_f/omega*2*pi,flow_RRF_3_nominal.pp001r5_fft,'DisplayName','RRF $p_{rms}$')
% loglog(flow_unst_3_nominal.pp001r5_fft_f/omega*2*pi,flow_unst_3_nominal.pp001r5_fftThick,'DisplayName','Unsteady Thick')
% loglog(flow_unst_3_nominal.pp001r5_fft_f/omega*2*pi,flow_unst_3_nominal.pp001r5_fftLoad,'DisplayName','Unsteady Loading')
% loglog(flow_RRF_3_nominal.pp001r5_fft_f/omega*2*pi,flow_RRF_3_nominal.pp001r5_fftThick,'DisplayName','RRF Thick')
% loglog(flow_RRF_3_nominal.pp001r5_fft_f/omega*2*pi,flow_RRF_3_nominal.pp001r5_fftLoad,'DisplayName','RRF Loading')

a = patch(pa1X/omega*2*pi,paY,'r','DisplayName','Inaudible frequencies');
b = patch(pa2X/omega*2*pi,paY,'r','HandleVisibility','off');
alpha(a,alpha_val);
alpha(b,alpha_val);
legend('Location','northeast')
title('')
xlabel('$f/f_{rot}$ [-]')
ylabel('p [Pa]')
xline(f/omega*2*pi,'-.',num2str(f/omega*2*pi),'DisplayName','Rotation frequency','LineWidth',2)
xline(Nb*f/omega*2*pi,'--',num2str(Nb*f/omega*2*pi),'DisplayName','Blade passing frequency','LineWidth',2)
xlim([0.7,700])
ylim([1e-6,0.2])
% subplot(2,1,2)
% loglog(flow_unst_3_nominal.CF_f/omega*2*pi,flow_unst_3_nominal.CF_fft,'DisplayName','CF')
% hold on;
% loglog(flow_unst_3_nominal.CL_f/omega*2*pi,flow_unst_3_nominal.CL_fft,'DisplayName','CL')
% legend
% xlim([1,10000])
exportgraphics(fig_fft,imagesPath+"fft_acoustic_pressure.emf")

% thickness and loading
fig_fft_thick_load = figure('Position',[100,100,pixel_x,pixel_y]);
loglog(flow_unst_3_nominal.pp001r5_fft_f/omega*2*pi,flow_unst_3_nominal.pp001r5_fftLoad,'b','DisplayName','Loading unsteady')
hold on;
loglog(flow_RRF_3_nominal.pp001r5_fft_f/omega*2*pi,flow_RRF_3_nominal.pp001r5_fftLoad,'r','DisplayName','Loading RRF')
loglog(flow_unst_3_nominal.pp001r5_fft_f/omega*2*pi,flow_unst_3_nominal.pp001r5_fftThick,'g','DisplayName','Thickness unsteady')
loglog(flow_RRF_3_nominal.pp001r5_fft_f/omega*2*pi,flow_RRF_3_nominal.pp001r5_fftThick,'m','DisplayName','Thickness RRF')
a = patch(pa1X/omega*2*pi,paY,'r','DisplayName','Inaudible frequencies');
b = patch(pa2X/omega*2*pi,paY,'r','HandleVisibility','off');
alpha(a,alpha_val);
alpha(b,alpha_val);
legend('Location','northeast')
title('')
xlabel('$f/f_{rot}$ [-]')
ylabel('p [Pa]')
xline(f/omega*2*pi,'-.',num2str(f/omega*2*pi),'DisplayName','Rotation frequency','LineWidth',2)
xline(Nb*f/omega*2*pi,'--',num2str(Nb*f/omega*2*pi),'DisplayName','Blade passing frequency','LineWidth',2)
xlim([0.7,700])
ylim([1e-6,0.2])
exportgraphics(fig_fft_thick_load,imagesPath+"fft_thick_load.emf")


% STROUHAL
fig_Strouhal = figure('Position',[100,100,pixel_x,pixel_y]);
loglog(flow_unst_3_nominal.pp001r5_PSD_St,flow_unst_3_nominal.pp001r5_PSD,'DisplayName','Nominal at 0')
hold on;
loglog(flow_RRF_3_nominal.pp001r5_PSD_St,flow_RRF_3_nominal.pp001r5_PSD,'DisplayName','Nominal at 90')
legend
title('')
xlabel('St [-]')
ylabel('p [Pa]')
xlim([20,450]*L/U)
xline(Nb*St,'--',num2str(Nb*St),'DisplayName','Nb/T','LineWidth',2)
% exportgraphics(fig_Strouhal,imagesPath+"Strouhal_acoustic_pressure.emf")
