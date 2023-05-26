%{
plot the dynamic CMz for unsteady simulations. 
Set the folders where you have the results in the simulationsFolderPath
string variable, and then set the name of the simulations in the
simulationsNames string variable
%}
%% init
clear; close all; clc; 
matlabCodesPath = "C:\Users\marco\OneDrive - Politecnico di Milano\MAGISTRALE\QuartoSemestre\Aeroacoustics\aeroacoustic project\matlabscripts";

% select simulations to extract
simulationsFolderPath = [   "E:\UNI - fisso\aeroacustica\unsteady_3_profili_deform6\";
                            "E:\UNI - fisso\aeroacustica\unsteady_3_profili_deform2\";
                            "E:\UNI - fisso\aeroacustica\unsteady_3_profili_nominal\"];

% "E:\UNI - fisso\aeroacustica\unsteady_3_profili_deform6\";
% "\\wsl.localhost\Ubuntu-20.04\home\marco\unsteady_3_profili_deform6\";
% "\\wsl.localhost\Ubuntu-20.04\home\marco\unsteady_3_profili_deform1\";
% "\\wsl.localhost\Ubuntu-20.04\home\marco\unsteady_1_profilo_deform2\";
% "\\wsl.localhost\Ubuntu-20.04\home\marco\unsteady_3_profili_deform7\"

% settare questo per i nomi nei plot
simulationsNames = ["3p deform6","3p deform2","3p nominal","1p deform7"]; 
matlab_graphics;

%% extract data
colors = ['r','b','k','m'];
omega = 21.33;
% check how to extract automatically the value for dt
dt_def = 8e-4;

for j = 1:length(simulationsFolderPath)
    if simulationsFolderPath(j) == "E:\UNI - fisso\aeroacustica\unsteady_3_profili_deform2\"
        dt(j) = 1e-4;
    else
        dt(j) = 2e-4;
    end

    history_files = dir(simulationsFolderPath(j)+"history*");
    len = 0;
    history = [];
    if length(history_files)~=1
        for i = 1:length(history_files)
    
            h_p = readmatrix(simulationsFolderPath(j)+history_files(i).name);
            if simulationsFolderPath(j)+history_files(i).name == "E:\UNI - fisso\aeroacustica\unsteady_3_profili_deform6\history_02334.dat"
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
            end
        end
    else
           
        history = readmatrix(simulationsFolderPath(j)+history_files.name);
        iterations_lost(j) = 0;
        time_loss(j) = 0;
    end
    %% adjust the coefficient   
    
    final_phase(j) = rad2deg(length(history(:,1))*dt(j)+iterations_lost(j)*dt(j))*omega;
    offset = 0; % offset phase in °, if 0 then it takes a single rotation from multiples of 120°
    iter_to_remove(j) = 0; %round(deg2rad(mod(final_phase(j)+offset,120))/omega/dt(j));
    
    history = history(1:end-iter_to_remove(j),:);

    % extract quantities
    t_vec{j} = [1:length(history(:,1))]*dt(j);
    angle{j} = rad2deg(t_vec{j} * omega);

    CMz{j} = history(:,11)/0.75/0.075;
    inner_iters{j} = history(:,1);
    rms_rho{j} = history(:,2);
    CD{j} = history(:,9);
    CL{j} = history(:,10);
    %% post process
    time_ratio = dt_def/dt(j);
    N_iter_per_round(j) = 368 * time_ratio;
    N_rounds_to_average(j) = floor(length(history(:,11))/N_iter_per_round(j));
    
    
    if N_rounds_to_average(j) ~= 0
        N_rounds_to_average(j) =1;
        starting_iter = size(history,1)-N_iter_per_round(j)*N_rounds_to_average(j);
        CMz_cycle_avg{j} = mean(CMz{j}(starting_iter + 1:starting_iter+N_iter_per_round(j)));
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

% RRF_CMz_trans_3 = 0.0609;
% RRF_CMz_trans_1 = 0.0947;
% RRF_CMz_mach = 0.147;

%% plots
figure('Position',[100,100,600,400])
subplot(2,1,1)
for j = 1:length(simulationsFolderPath)
    plot(t_vec{j}+time_loss(j),CMz{j},[colors(j),'-'],'DisplayName',"Unsteady "+simulationsNames(j))
    hold on; 
    % yline(RRF_CMz_trans_1,'k--',"RRF trans ="+num2str(RRF_CMz_trans_1),'HandleVisibility','off')
    if ~isempty(CMz_cycle_avg{j})
        yline(CMz_cycle_avg{j}(end),[colors(j),'--'],"CMz avg = "+num2str(CMz_cycle_avg{j}(end)),'HandleVisibility','off')
        xline(t_vec{j}(end-N_rounds_to_average(j)*N_iter_per_round(j)),[colors(j),'-'],"AFH",'HandleVisibility','off')
    end
end
legend
subplot(2,1,2)
for j = 1:length(simulationsFolderPath)
    plot(t_vec{j}+time_loss(j),rms_rho{j},[colors(j),'-'],'DisplayName',"Unsteady "+simulationsNames(j))
    hold on; 
end
yline(-8,'k--',"target",'HandleVisibility','off')
legend

