%{
plot the dynamic CMz for unsteady simulations. 
Set the folders where you have the results in the simulationsFolderPath
string variable, and then set the name of the simulations in the
simulationsNames string variable
%}
%% init
clear; close all; clc; 
user = 'marco';
user_settings;

%% select simulations to extract
%unsteady
unstSimulationsFolderPath = ["\\wsl.localhost\Ubuntu-20.04\home\marco\unsteady_3_profili_surface_1p\";
                            "E:\UNI - fisso\aeroacustica\unsteady_1_profilo_nominal_dt2\";
                            "E:\UNI - fisso\aeroacustica\unsteady_3_profili_nominal\";];

simulationsNames = ["Hyb nom","1p Nom","3p Nom", "3p nom"]; 


% rrf mach
rrfMachSimulationsFolderPath = ["C:\Users\marco\OneDrive - Politecnico di Milano\MAGISTRALE\QuartoSemestre\Aeroacoustics\aeroacoustic project\Simulazioni\RRF\RRF_1P_MACH_0\";
                         "C:\Users\marco\OneDrive - Politecnico di Milano\MAGISTRALE\QuartoSemestre\Aeroacoustics\aeroacoustic project\Simulazioni\RRF\RRF_1P_MACH_90\";
                          "C:\Users\marco\OneDrive - Politecnico di Milano\MAGISTRALE\QuartoSemestre\Aeroacoustics\aeroacoustic project\Simulazioni\RRF\RRF_1P_MACH_180\";
                         "C:\Users\marco\OneDrive - Politecnico di Milano\MAGISTRALE\QuartoSemestre\Aeroacoustics\aeroacoustic project\Simulazioni\RRF\RRF_1P_MACH_270\"];


rrfMachNames = ["0", "90", "180", "270"];


% rrf trans
rrfTransSimulationsFolderPath = ["C:\Users\marco\OneDrive - Politecnico di Milano\MAGISTRALE\QuartoSemestre\Aeroacoustics\aeroacoustic project\Simulazioni\RRF\RRF_1P_TRANS_0\";
                             "C:\Users\marco\OneDrive - Politecnico di Milano\MAGISTRALE\QuartoSemestre\Aeroacoustics\aeroacoustic project\Simulazioni\RRF\RRF_1P_TRANS_90\";
                             "C:\Users\marco\OneDrive - Politecnico di Milano\MAGISTRALE\QuartoSemestre\Aeroacoustics\aeroacoustic project\Simulazioni\RRF\RRF_1P_TRANS_180\";
                             "C:\Users\marco\OneDrive - Politecnico di Milano\MAGISTRALE\QuartoSemestre\Aeroacoustics\aeroacoustic project\Simulazioni\RRF\RRF_1P_TRANS_270\";];
rrfTransNames = ["0", "90", "180", "270"];

% "E:\UNI - fisso\aeroacustica\unsteady_1_profilo_deform2\";
% "E:\UNI - fisso\aeroacustica\unsteady_3_profili_deform6\";
% "\\wsl.localhost\Ubuntu-20.04\home\marco\unsteady_3_profili_deform6\";
% "\\wsl.localhost\Ubuntu-20.04\home\marco\unsteady_3_profili_deform1\";
% "\\wsl.localhost\Ubuntu-20.04\home\marco\unsteady_1_profilo_deform2\";
% "\\wsl.localhost\Ubuntu-20.04\home\marco\unsteady_3_profili_deform7\"

% settare questo per i nomi nei plot
matlab_graphics;

%% extract data - unsteady
colors = ['r','b','k','m'];
omega = 21.33;
% check how to extract automatically the value for dt
dt_def = 8e-4;

for j = 1:length(unstSimulationsFolderPath)
    if unstSimulationsFolderPath(j) == "E:\UNI - fisso\aeroacustica\unsteady_3_profili_deform2\" || unstSimulationsFolderPath(j) =="C:\Users\marco\OneDrive - Politecnico di Milano\MAGISTRALE\QuartoSemestre\Aeroacoustics\aeroacoustic project\Simulazioni\Unsteady\Unsteady_1_profilo_nominal_dt1\"
        dt(j) = 1e-4;
    elseif unstSimulationsFolderPath(j) == "E:\UNI - fisso\aeroacustica\unsteady_1_profilo_nominal\" || unstSimulationsFolderPath(j) =="C:\Users\marco\OneDrive - Politecnico di Milano\MAGISTRALE\QuartoSemestre\Aeroacoustics\aeroacoustic project\Simulazioni\Unsteady\Unsteady_1_profilo_nominal_dt4\"
        dt(j) = 4e-4;
    else
        dt(j) = 2e-4;
    end

    history_files = dir(unstSimulationsFolderPath(j)+"history*");
    len = 0;
    history = [];
    if length(history_files)~=1
        for i = 1:length(history_files)
    
            h_p = readmatrix(unstSimulationsFolderPath(j)+history_files(i).name);
            if unstSimulationsFolderPath(j)+history_files(i).name == "E:\UNI - fisso\aeroacustica\unsteady_3_profili_deform6\history_02334.dat"
                h_p = h_p(1:2:end,:);
            end
            restart_iter = str2double(erase(history_files(i).name,["history_",".dat"]));
            if restart_iter < length(history)
                iter_to_overwrite = size(history,1)-restart_iter;
            else
                iter_to_overwrite = 0;
            end
            history = history(1:end-iter_to_overwrite,:);
            history(len+1-iter_to_overwrite:len+length(h_p)-iter_to_overwrite,:) = h_p;
            
            len = length(history);
            

            if i == length(history_files)

                iterations_lost(j) = str2double(erase(history_files(end).name,["history_",".dat"]))+length(h_p)-len;
                if iterations_lost(j)  > 0
                    time_loss(j) = (iterations_lost(j)-368*dt_def/dt(j)) * dt(j) ;  
                else
                    time_loss(j) = 0;
                end
                angle_loss(j) = time_loss(j)*omega;
            end
        end
    else
           
        history = readmatrix(unstSimulationsFolderPath(j)+history_files.name);
        if contains(history_files.name,'_')
            iterations_lost(j) = str2double(erase(history_files(end).name,["history_",".dat"]));             
            time_loss(j) = (iterations_lost(j)-1472) * dt(j);
          %  time_loss(j) = (iterations_lost(j)) * dt(j);
            angle_loss(j) = rad2deg(time_loss(j)*omega);
        else
            iterations_lost(j) = 0;
            time_loss(j) = 0;
            angle_loss(j) = 0;
        end
        
        
    end
    %% adjust the coefficient   
    
    final_phase(j) = rad2deg(length(history(:,1))*dt(j)+iterations_lost(j)*dt(j))*omega;
    offset = 0; % offset phase in °, if 0 then it takes a single rotation from multiples of 120°
    iter_to_remove(j) = 0; %round(deg2rad(mod(final_phase(j)+offset,120))/omega/dt(j));
    
    history = history(1:end-iter_to_remove(j),:);

    % extract quantities
    t_vec{j} = [1:length(history(:,1))]*dt(j);
    angle{j} = rad2deg(t_vec{j} * omega);

    if unstSimulationsFolderPath(j) == "E:\UNI - fisso\aeroacustica\unsteady_1_profilo_nominal\" || unstSimulationsFolderPath(j) == "C:\Users\marco\OneDrive - Politecnico di Milano\MAGISTRALE\QuartoSemestre\Aeroacoustics\aeroacoustic project\Simulazioni\Unsteady\Unsteady_1_profilo_nominal_dt4\"
        CMz{j} = history(:,11);
        CD{j} = history(:,9);
        CL{j} = history(:,10);
    else
        CMz{j} = history(:,11)/0.75/0.075;
        CD{j} = history(:,9)/0.075;
        CL{j} = history(:,10)/0.075;
    end
    inner_iters{j} = history(:,1);
    rms_rho{j} = history(:,2);
    
    
    %% post process
    time_ratio = dt_def/dt(j);
    N_iter_per_round(j) = 368 * time_ratio;
    N_rounds_to_average(j) = floor(length(history(:,11))/N_iter_per_round(j));
    
    
    if N_rounds_to_average(j) ~= 0
        N_rounds_to_average(j) =1;
        starting_iter = size(history,1)-N_iter_per_round(j)*N_rounds_to_average(j);
        CMz_cycle_avg{j} = mean(CMz{j}(starting_iter + 1:starting_iter+N_iter_per_round(j)));
        CD_cycle_avg{j} = mean(CD{j}(starting_iter + 1:starting_iter+N_iter_per_round(j)));
        CL_cycle_avg{j} = mean(CL{j}(starting_iter + 1:starting_iter+N_iter_per_round(j)));
        CMz_excursion{j} = max(CMz{j}(starting_iter + 1:starting_iter+N_iter_per_round(j)))-min(CMz{j}(starting_iter + 1:starting_iter+N_iter_per_round(j)));
        CMz_rms{j} = rms(CMz{j}(starting_iter + 1:starting_iter+N_iter_per_round(j))-mean(CMz{j}(starting_iter + 1:starting_iter+N_iter_per_round(j))));
    
        fprintf('CMz for %s = %.6d\n',simulationsNames(j),CMz_cycle_avg{j})
        fprintf('CMz excursion in last cycle for %s = %.3d\n',simulationsNames(j),CMz_excursion{j})
        fprintf('CMz rms in last cycle for %s = %.3d\n\n',simulationsNames(j),CMz_rms{j})
    else
        CMz_cycle_avg{j} = [];
        CMz_excursion{j} = [];
        CMz_rms{j} = [];
    end

    % remove variables so they don't overwrite badly
    clearvars history h_p
end

%% extract data - RRF MACH
RRF_CMz_mach = [];
for j = 1:length(rrfMachSimulationsFolderPath)
RRF_mach = readmatrix(rrfMachSimulationsFolderPath(j)+"history.dat");
RRF_CMz_mach(j) = RRF_mach(end,end-1);
RRF_CL_mach(j) = RRF_mach(end,end-2);
RRF_CD_mach(j) = RRF_mach(end,end-3);
end

%% extract data - RRF TRANS
RRF_CMz_trans = [];
for j = 1:length(rrfTransSimulationsFolderPath)
RRF_trans = readmatrix(rrfTransSimulationsFolderPath(j)+"history.dat");
RRF_CMz_trans(j) = RRF_trans(end,end-1);
RRF_CL_trans(j) = RRF_trans(end,end-2);
RRF_CD_trans(j) = RRF_trans(end,end-3);
end

%% plots
markers = 'sd^o';
angRRF = [0, 90, 180, 270];
round = 2;

% CMz
fig_CMz_RRF_vs_unsteady =figure('Position',[100,100,900,600]);
for j = 1:length(unstSimulationsFolderPath)
    plot(angle{j}+angle_loss(j)-(round-1)*360,CMz{j},[colors(j),'-'],'DisplayName',"Unsteady "+simulationsNames(j))
    hold on; 
%     yline(RRF_CMz_trans_1,'k--',"RRF trans ="+num2str(RRF_CMz_trans_1),'HandleVisibility','off')
    if ~isempty(CMz_cycle_avg{j})
        yline(CMz_cycle_avg{j}(end),[colors(j),'--'],'DisplayName',"CMz avg = "+num2str(CMz_cycle_avg{j}(end)))
%         xline(t_vec{j}(end-N_rounds_to_average(j)*N_iter_per_round(j)),[colors(j),'-'],"AFH",'HandleVisibility','off')
    end
end
% if ~isempty(RRF_CMz_mach)
%     scatter(angRRF,RRF_CMz_mach,'gs',"DisplayName","RRF Mach")  
% end
% if ~isempty(RRF_CMz_trans)
%     scatter(angRRF,RRF_CMz_trans,'b^',"DisplayName","RRF Trans")    
% end
legend('Location','northeast')
ylim([-0.3, 0.5])
xlim([0, 1080])
xlabel('Angle [°]')
ylabel('CMz [-]')
% exportgraphics(fig_CMz_RRF_vs_unsteady,imagesPath+"22_CMz_nominal_vs_deform6.emf")


%% CD
fig_CD_RRF_vs_unsteady =figure('Position',[100,100,900,600]);
for j = 1:length(unstSimulationsFolderPath)
    plot(angle{j}+angle_loss(j)-(round-1)*360,CD{j},[colors(j),'-'],'DisplayName',"Unsteady "+simulationsNames(j))
    hold on; 
    % yline(RRF_CMz_trans_1,'k--',"RRF trans ="+num2str(RRF_CMz_trans_1),'HandleVisibility','off')
    if ~isempty(CD_cycle_avg{j})
        yline(CD_cycle_avg{j}(end),[colors(j),'--'],'DisplayName',"CD avg = "+num2str(CD_cycle_avg{j}(end)))
%         xline(t_vec{j}(end-N_rounds_to_average(j)*N_iter_per_round(j)),[colors(j),'-'],"AFH",'HandleVisibility','off')
    end
end
if ~isempty(RRF_CD_mach)
    scatter(angRRF,RRF_CD_mach,'gs',"DisplayName","RRF Mach")
end
if ~isempty(RRF_CD_trans)
    scatter(angRRF,RRF_CD_trans,'b^',"DisplayName","RRF Trans")
end
legend('Location','northeast')
% ylim([-0.3, 1])
xlim([0, 1080])
xlabel('Angle [°]')
ylabel('CD [-]')
% exportgraphics(fig_CD_RRF_vs_unsteady,imagesPath+"11_CD_RRF_vs_unsteady.emf")


%CL
fig_CL_RRF_vs_unsteady =figure('Position',[100,100,900,600]);
for j = 1:length(unstSimulationsFolderPath)
    plot(angle{j}+angle_loss(j)-(round-1)*360,CL{j},[colors(j),'-'],'DisplayName',"Unsteady "+simulationsNames(j))
    hold on; 
    % yline(RRF_CMz_trans_1,'k--',"RRF trans ="+num2str(RRF_CMz_trans_1),'HandleVisibility','off')
    if ~isempty(CL_cycle_avg{j})
        yline(CL_cycle_avg{j}(end),[colors(j),'--'],'DisplayName',"CL avg = "+num2str(CL_cycle_avg{j}(end)))
%         xline(t_vec{j}(end-N_rounds_to_average(j)*N_iter_per_round(j)),[colors(j),'-'],"AFH",'HandleVisibility','off')
    end
end
if ~isempty(RRF_CL_mach)
    scatter(angRRF,RRF_CL_mach,'gs',"DisplayName","RRF Mach")
end
if ~isempty(RRF_CL_trans)
    scatter(angRRF,RRF_CL_trans,'b^',"DisplayName","RRF Trans")
end
legend
% ylim([-0.3, 1])
xlim([0, 1080])
xlabel('Angle [°]')
ylabel('CL [-]')
% exportgraphics(fig_CMz_RRF_vs_unsteady,imagesPath+"11_CL_RRF_vs_unsteady.emf")

% % timestep convergence
% idx_cut = [4520/4,4520/2,4520];
% 
% fig_timeStepConvergence = figure;%('Position',[100,100,600,400]);
% for j = 1:length(unstSimulationsFolderPath)
%     plot(t_vec{j}(idx_cut(j)+1:idx_cut(j)+N_iter_per_round(j))+time_loss(j),CMz{j}(idx_cut(j)+1:idx_cut(j)+N_iter_per_round(j)),[colors(j),'-'],'DisplayName',"Unsteady "+simulationsNames(j))
%     hold on; grid on;
% end
% xlim([min(t_vec{1}(idx_cut(1)+1:idx_cut(1)+N_iter_per_round(1))),max(t_vec{1}(idx_cut(1)+1:idx_cut(1)+N_iter_per_round(1)))])
% xlabel('Time [s]')
% ylabel('CMz [-]')
% legend
% exportgraphics(fig_timeStepConvergence,imagesPath+"08_timeStepConv.emf")

