%{

Observer location generator

%}
clear; close all; clc; 
folder_wsl = "E:\UNI - fisso\aeroacustica\unsteady_3_profili_deform1\observer_05\";
folder_drive = "C:\Users\marco\OneDrive - Politecnico di Milano\MAGISTRALE\QuartoSemestre\Aeroacoustics\aeroacoustic project\varie\";
user = 'marco';
user_settings;
list = dir(folder_wsl+"pp_FWH_*");
for i = 1:length(list)
    pp{i} = readmatrix(folder_wsl+list(i).name);
end

% extract locations
obs_loc = readmatrix(folder_wsl+"Observer_Locations.dat");
profile = readmatrix(folder_drive+"NACA0021.txt");
profile1 = profile([2:801], :);
profile2 = profile([803:1602], :);
profile3 = profile([1604:end], :);


%% plots
% pressure fluctuations
% figure
% hold on;
% for i = 1:length(pp)
% plot(pp{i}(:,1),pp{i}(:,2))
% end


% observer locations
fig_observer_loc = figure('Position',[100,100,600,500]);
scatter(obs_loc(:,1),obs_loc(:,2),'DisplayName','Observer locations')
hold on; grid minor;
patch(profile1(:,1),profile1(:,2),'black','DisplayName','Airfoil')
patch(profile2(:,1),profile2(:,2),'black','HandleVisibility','off')
patch(profile3(:,1),profile3(:,2),'black','HandleVisibility','off')
% plot(newPoints(:,1),newPoints(:,2))
axis equal
legend
exportgraphics(fig_observer_loc,folder_drive+"observer_location.png")


%% generate observers

N = 120;
r = [5, 10, 20];
th = linspace(0,360,N+1);
th = th(1:end-1);

for i = 1:length(r)
    newPoints{i} = r(i)*[cosd(th);sind(th);zeros(1,N)]';
    fid = fopen(['Observer_Locations_r',num2str(r(i)),'_',num2str(N),'.dat'],'w');
    fprintf(fid,'%d\n',N);
    for j = 1: size(newPoints{i},1)
        fprintf(fid,'%.4f %.4f %d\n', newPoints{i}(j,1),newPoints{i}(j,2),newPoints{i}(j,3))
    end
    fclose(fid);
end
% observer locations
fig_observer_loc = figure('Position',[100,100,600,500]);
hold on; grid on;
scatter(newPoints{1}(1:4:end,1),newPoints{1}(1:4:end,2),'DisplayName','New observer locations 5')
scatter(newPoints{2}(1:2:end,1),newPoints{2}(1:2:end,2),'DisplayName','New observer locations 10')
scatter(newPoints{3}(:,1),newPoints{3}(:,2),'DisplayName','New observer locations 20')
patch(profile1(:,1),profile1(:,2),'black','DisplayName','Airfoil')
patch(profile2(:,1),profile2(:,2),'black','HandleVisibility','off')
patch(profile3(:,1),profile3(:,2),'black','HandleVisibility','off')
% plot(newPoints(:,1),newPoints(:,2))
axis equal
legend
exportgraphics(fig_observer_loc,imagesPath+"10_timeStepConv.emf")