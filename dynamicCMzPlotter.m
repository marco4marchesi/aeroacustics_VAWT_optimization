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

%"\\wsl.localhost\Ubuntu-20.04\home\marco\unsteady_3_profili_deform6\";
%"\\wsl.localhost\Ubuntu-20.04\home\marco\unsteady_3_profili_deform1\";
%"\\wsl.localhost\Ubuntu-20.04\home\marco\unsteady_1_profilo_deform2\";
%"\\wsl.localhost\Ubuntu-20.04\home\marco\unsteady_3_profili_drag_fixedMoment\"
%"3p drag fixedMoment",

% settare questo per i nomi nei plot
simulationsNames = ["3p deform6","3p drag fixM","3p deform1","3p nominal"]; 
matlab_graphics;

%% extract data
colors = ['r','b','k','m'];
omega = 21.33;
% check how to extract automatically the value for dt
dt_def = 8e-4;
dt = 2e-4;
for j = 1:length(simulationsFolderPath)
    
    history_files = dir(simulationsFolderPath(j)+"history*");
    len = 0;
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
            
            history(len+1:len+length(h_p),:) = h_p;
            len = length(history);
            iterations_lost(j) = str2double(erase(history_files(end).name,["history_",".dat"]))+length(h_p)-len;
            time_loss(j) = (iterations_lost(j)-1472) * dt ;  
            
        end
    else
           
        history = readmatrix(simulationsFolderPath(j)+history_files.name);
        iterations_lost(j) = 0;
        time_loss(j) = 0;
    end
    %% adjust the coefficient   
    
    final_phase(j) = rad2deg(length(history(:,1))*dt+iterations_lost(j)*dt)*omega;
    offset = 50; % offset phase in °, if 0 then it takes a single rotation from multiples of 120°
    iter_to_remove(j) = round(deg2rad(mod(final_phase(j)+offset,120))/omega/dt);
    
    history = history(1:end-iter_to_remove(j),:);

    % extract quantities
    t_vec{j} = [1:length(history(:,1))]*dt;
    angle{j} = rad2deg(t_vec{j} * omega);

    CMz{j} = history(:,11)/0.75/0.075;
    inner_iters{j} = history(:,1);
    rms_rho{j} = history(:,2);
    CD{j} = history(:,9);
    CL{j} = history(:,10);
    %% post process
    time_ratio = dt_def/dt;
    N_iter_per_round = 368 * time_ratio;
    N_rounds_to_average(j) = floor(length(history(:,11))/N_iter_per_round);
    
    
    if N_rounds_to_average(j) ~= 0
        N_rounds_to_average(j) =1;
        starting_iter = size(history,1)-N_iter_per_round*N_rounds_to_average(j);
        for i =  1: N_rounds_to_average(j)
            CMz_cycle_avg{j}(i) = mean(CMz{j}((i-1)*N_iter_per_round+starting_iter + 1:i*N_iter_per_round+starting_iter));
        end
    else
        CMz_cycle_avg{j} = [];
    end

    % remove variables so they don't overwrite badly
    clearvars history h_p
end

% RRF_CMz_trans_3 = 0.0609;
% RRF_CMz_trans_1 = 0.0947;
% RRF_CMz_mach = 0.147;

%% plots
figure('Position',[100,100,600,400])
for j = 1:length(simulationsFolderPath)
    plot(t_vec{j}+time_loss(j),CMz{j},[colors(j),'-'],'DisplayName',"Unsteady "+simulationsNames(j))
    hold on; 
    % yline(RRF_CMz_trans_1,'k--',"RRF trans ="+num2str(RRF_CMz_trans_1),'HandleVisibility','off')
    if ~isempty(CMz_cycle_avg{j})
        yline(CMz_cycle_avg{j}(end),[colors(j),'--'],"CMz avg = "+num2str(CMz_cycle_avg{j}(end)),'HandleVisibility','off')
        xline(t_vec{j}(end-N_rounds_to_average(j)*N_iter_per_round),[colors(j),'-'],"AFH",'HandleVisibility','off')
    end
end
legend

% 
% figure('Position',[100,100,600,400])
% for j = 3
%     plot(angle{j},CL{j},'r-','DisplayName',"CL "+simulationsNames(j))
%     hold on; 
%     plot(angle{j},CD{j},'b-','DisplayName',"CD "+simulationsNames(j))
%     plot(angle{j},CMz{j},'k-','DisplayName',"CMz "+simulationsNames(j))
% end
% xline(360,'--','DisplayName','360°')
% xline(450,'--','DisplayName','450°')
% grid on;
% xlabel('Angle [°]')
% ylabel('Coeff [-]')
% legend
