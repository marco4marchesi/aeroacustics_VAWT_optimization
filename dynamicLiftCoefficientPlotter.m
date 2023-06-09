%{

Matlab code for extracting data from history files of multiple angles simulations.

--------------------------------------------------------------------------
Author: Marco Marchesi
--------------------------------------------------------------------------

+   this code aims to extract the CD, CL, and CM ---curves--- w.r.t. the angle of 
    attack from the simulations and to automatically save them in a vector that 
    contains the number of elements in the mesh

+   this is done for checking the convergence of the meshes, we need to
    investigate whether the mesh returns a better result w.r.t. the precedent
    tried mesh

---------------------- HOW TO USE THIS CODE: -----------------------------

if it's the first time you use this code follow the following steps before
running the code, otherwise it will return errors:

+   set the "matlabCodesPath" folder equal to the path where you have this
    script. 
    To do so you have to change the line in the INIT section  (lines 20 -
    40) where you see your name. If you don't see your name just copy-paste
    one "if" and set your "matlabCodesPath" variable.

+   set the "simulationsFolderPath" to the folder where you have all you
    simulations saved. The last folder should be "Simulation/". if it's not
    like that try asking for help to the author.


%}

%% select user

% user: set who is running the code so that the folder is chosen:
user = "doppio fisso"; % choices: "doppio fisso" "luca" ...

if user == "doppio fisso"
    matlabCodesPath = "C:\Users\marco\OneDrive - Politecnico di Milano\MAGISTRALE\QuartoSemestre\Aeroacoustics\aeroacoustic project\matlabscripts";
    simulationsFolderPath = "\\wsl.localhost\Ubuntu-20.04\home\marco\test_script";
%     simulationsFolderPath = "C:\Users\marco\Desktop\UNI\2 MAGISTRALE\CFD\CFD PROJECT\progetto_CFD\Simulations\"; % locale
end

if user == "doppio portatile"
    matlabCodesPath = "C:\Users\marco\Desktop\tutto\UNI\2 MAGISTRALE\CFD\CFD PROJECT\progetto_CFD\Codes\matlabCodes";
    simulationsFolderPath = "C:\Users\marco\OneDrive - Politecnico di Milano\MAGISTRALE\TerzoSemestre\CFD\PROGETTO CFD DRIVE\SIMULATIONS DRIVE\Simulations\";
%     simulationsFolderPath = "C:\Users\marco\Desktop\UNI\2 MAGISTRALE\CFD\CFD PROJECT\progetto_CFD\Simulations\"; % locale
end

if user == "luca"
    matlabCodesPath = "C:/Users/lucag/Desktop/Universita/Magistrale Secondo Anno/Computational_fluid_dynamics/Progetto_CFD/progetto_CFD/Codes/matlabCodes\";
    simulationsFolderPath = "C:/Users/lucag/Desktop/Universita/Magistrale Secondo Anno/Computational_fluid_dynamics/Progetto_CFD/progetto_CFD/Simulations\";
end

%% init 
clearvars -except matlabCodesPath simulationsFolderPath; 
close all; clc;

% add matlab functions to the path
% rmpath(matlabCodesPath+"/convergence_analysis")
% rmpath(matlabCodesPath+"/polar_plotter")
% rmpath(matlabCodesPath+"/farfield_analysis")
% addpath(matlabCodesPath)
% addpath(matlabCodesPath+"/dynamic_lift_coefficient")
% addpath(matlabCodesPath+"/utilitiesFunctions")
% move to simulation folder
cd(simulationsFolderPath)
% matlab_graphics;




%% ------------------------------------ CHOSE SIMULATION (FOLDER) -------------------------------------- %%
mainFolder = "";%["TC1";"TC2";"TC1_LM";"TC2_LM";"TC3_LM";"TC4_LM"];
simuFolder = '"";%/HISTORY'; % if some errors occur try using single apices

% do you want to plot ONLY THE LAST PERIOD?
casePeriod = 'first'; % options: {'last', 'first', ''}



%% extract Zanotti dataset


% extract_ZANO_CL_exp = readtable(simulationsFolderPath+"../plotDigitizer/DynamicCL_Zanotti_experimental.csv","Delimiter",';');
% extract_ZANO_CL_CFD_720TS = readtable(simulationsFolderPath+"../plotDigitizer/DynamicCL_Zanotti_CFD_G1_720TS.csv","Delimiter",';');
% extract_ZANO_CM_exp = readtable(simulationsFolderPath+"../plotDigitizer/DynamicCM_Zanotti_experimental.csv","Delimiter",';');
% extract_ZANO_CM_CFD_720TS = readtable(simulationsFolderPath+"../plotDigitizer/DynamicCM_Zanotti_CFD_G1_720TS.csv","Delimiter",';');


% for i = 1: size(extract_ZANO_CL_CFD_720TS,1)
%     dynamicCL_Zanotti_CFD.alpha(i,1) = str2double(replace(extract_ZANO_CL_CFD_720TS{i,1}{1},",",".")) ;
%     dynamicCL_Zanotti_CFD.CL(i,1) = str2double(replace(extract_ZANO_CL_CFD_720TS{i,2}{1},",","."));
% end


% extract area integral
% [~,idx_minzanocfd]  = min(dynamicCL_Zanotti_CFD.alpha);
% [~,idx_maxzanocfd]  = max(dynamicCL_Zanotti_CFD.alpha);

% if idx_minzanocfd<idx_maxzanocfd
%     area_upCFD = trapz(dynamicCL_Zanotti_CFD.alpha(idx_minzanocfd:idx_maxzanocfd),dynamicCL_Zanotti_CFD.CL(idx_minzanocfd:idx_maxzanocfd));
%     area_downCFD = trapz(dynamicCL_Zanotti_CFD.alpha(idx_maxzanocfd:end),dynamicCL_Zanotti_CFD.CL(idx_maxzanocfd:end));
% else
%     area_downCFD = trapz(dynamicCL_Zanotti_CFD.alpha(idx_maxzanocfd:idx_minzanocfd),dynamicCL_Zanotti_CFD.CL(idx_maxzanocfd:idx_minzanocfd));
%     area_upCFD = trapz(dynamicCL_Zanotti_CFD.alpha(1:idx_maxzanocfd),dynamicCL_Zanotti_CFD.CL(1:idx_maxzanocfd));
% end
% area_zanoCFD = area_upCFD-area_downCFD;




% % experiment
% for i = 1: size(extract_ZANO_CL_exp,1)
%     dynamicCL_Zanotti_exp.alpha(i,1) = str2double(replace(extract_ZANO_CL_exp{i,1}{1},",",".")) ;
%     dynamicCL_Zanotti_exp.CL(i,1) = str2double(replace(extract_ZANO_CL_exp{i,2}{1},",","."));
% end
% 
% % extract area integral
% [~,idx_minzanoexp]  = min(dynamicCL_Zanotti_exp.alpha);
% [~,idx_maxzanoexp]  = max(dynamicCL_Zanotti_exp.alpha); 
% 
% if idx_minzanoexp<idx_maxzanoexp
%     area_upEXP = trapz(dynamicCL_Zanotti_exp.alpha(idx_minzanoexp:idx_maxzanoexp),dynamicCL_Zanotti_exp.CL(idx_minzanoexp:idx_maxzanoexp));
%     area_downEXP = trapz(dynamicCL_Zanotti_exp.alpha(idx_maxzanoexp:end),dynamicCL_Zanotti_exp.CL(idx_maxzanoexp:end))+trapz(dynamicCL_Zanotti_exp.alpha(1:idx_minzanoexp),dynamicCL_Zanotti_exp.CL(1:idx_minzanoexp));
% else
%     area_downEXP = trapz(dynamicCL_Zanotti_exp.alpha(idx_maxzanoexp:idx_minzanoexp),dynamicCL_Zanotti_exp.CL(idx_maxzanoexp:idx_minzanoexp));
%     area_upEXP = trapz(dynamicCL_Zanotti_exp.alpha(1:idx_maxzanoexp),dynamicCL_Zanotti_exp.CL(1:idx_maxzanoexp));
% end
% area_zanoEXP = area_upEXP-area_downEXP;

% % cfd
% for i = 1: size(extract_ZANO_CM_CFD_720TS,1)
%     dynamicCM_Zanotti_CFD.alpha(i,1) = str2double(replace(extract_ZANO_CM_CFD_720TS{i,1}{1},",",".")) ;
%     dynamicCM_Zanotti_CFD.CM(i,1) = str2double(replace(extract_ZANO_CM_CFD_720TS{i,2}{1},",","."));
% end
% 
% % experiment
% for i = 1: size(extract_ZANO_CM_exp,1)
%     dynamicCM_Zanotti_exp.alpha(i,1) = str2double(replace(extract_ZANO_CM_exp{i,1}{1},",",".")) ;
%     dynamicCM_Zanotti_exp.CM(i,1) = str2double(replace(extract_ZANO_CM_exp{i,2}{1},",","."));
% end

%% cfd data
for kk = 1:size(mainFolder,1)
testcase = strcat(mainFolder(kk,:),simuFolder);
cd(testcase)


%% data extraction and post processing

% data extraction
listing = dir("*.csv");

history = struct('Inner_Iter',[],'CL',[],'CMz',[],'cauchyCL',[]);
for j = 1:length(listing)
currentHistory = csvDataLogExtractor(listing(j).name);
history.Inner_Iter = [history.Inner_Iter; currentHistory.Inner_Iter(currentHistory.Inner_Iter~=0)];
history.CL = [history.CL; currentHistory.CL(currentHistory.Inner_Iter~=0)];
history.CMz = [history.CMz; currentHistory.CMz(currentHistory.Inner_Iter~=0)];
history.cauchyCL = [history.cauchyCL; currentHistory.Cauchy_CL_(currentHistory.Inner_Iter~=0)];
end
firstHistoryStep = str2double(listing(1).name(end-8:end-4));

% simulation parameters (build alpha function)

A = 10;                             % [°] amplitude
omega = 13.5;                       % [rad/s] pitching pulsation of the simulation
alpha_mean = 10;                    % [°] mean angle of the simulation
chord = 1.00898;                    % [m] airfoil chord length
fsV = 68.05093;                     % [m/s] freestream velocity

switch mainFolder(kk,:) 
    case 'TC1'
        time_step(kk)= 0.02 * chord/fsV;        % [s] timestep of the simulation
        T_period(kk) = 1571;
    case 'TC2'
        time_step(kk)= 5*0.02 * chord/fsV;      % [s] timestep of the simulation
        T_period(kk) = 314;
    case 'TC1_LM'
        time_step(kk)= 0.02 * chord/fsV;        % [s] timestep of the simulation
        T_period(kk) = 1571;
    case 'TC2_LM'
        time_step(kk)= 5*0.02 * chord/fsV;      % [s] timestep of the simulation
        T_period(kk) = 314;
    case 'TC3_LM'
        time_step(kk)= 25*0.02 * chord/fsV;      % [s] timestep of the simulation
        T_period(kk) = 63;
    case 'TC4_LM'
        time_step(kk)= 5*0.02 * chord/fsV;      % [s] timestep of the simulation
        T_period(kk) = 314;
end



time_iter{:,kk} = 1:length(history.Inner_Iter);
inner_iter{:,kk} = history.Inner_Iter;
t{:,kk} = time_iter{:,kk} * time_step(kk)+firstHistoryStep*time_step(kk);
CL{:,kk} = history.CL;
CMz{:,kk} = history.CMz;
% compute alpha
alpha{:,kk} = A * sin(omega*(t{:,kk})) + alpha_mean; % [°] variable angle of the simulation
cauchyCL{:,kk} = history.cauchyCL;



% if you want to skip the first quarter
time_iter_offset = floor(T_period(kk)/4)+1; % to take the period from the maximum angle




    switch casePeriod
        case 'last'
    time_iter_plot{:,kk} = time_iter{:,kk}(end-T_period(kk):end);
    inner_iter_plot{:,kk} = inner_iter{:,kk}(end-T_period(kk):end);
    t_plot{:,kk} = t{:,kk}(end-T_period(kk):end);
    alpha_plot{:,kk} = alpha{:,kk}(end-T_period(kk):end);
    cauchyCL_plot{:,kk} = cauchyCL{:,kk}(end-T_period(kk):end);
    CL_plot{:,kk} = CL{:,kk}(end-T_period(kk):end);
    CMz_plot{:,kk} = CMz{:,kk}(end-T_period(kk):end);

    %     correction{:,kk} = correction{:,kk}(end-T_period(kk):end);

        case 'first'
    time_iter_plot{:,kk} = time_iter{:,kk}(time_iter_offset:T_period(kk)+time_iter_offset);
    inner_iter_plot{:,kk} = inner_iter{:,kk}(time_iter_offset:T_period(kk)+time_iter_offset);
    t_plot{:,kk} = t{:,kk}(time_iter_offset:T_period(kk)+time_iter_offset);
    alpha_plot{:,kk} = alpha{:,kk}(time_iter_offset:T_period(kk)+time_iter_offset);
    cauchyCL_plot{:,kk} = cauchyCL{:,kk}(time_iter_offset:T_period(kk)+time_iter_offset);
    CL_plot{:,kk} = CL{:,kk}(time_iter_offset-1:T_period(kk)+time_iter_offset-1);
    CMz_plot{:,kk} = CMz{:,kk}(time_iter_offset:T_period(kk)+time_iter_offset);
%     correction{:,kk} = correction{:,kk}(1:T_period(kk));
        case ''
    time_iter_plot{:,kk} = time_iter{:,kk};
    inner_iter_plot{:,kk} = inner_iter{:,kk};
    t_plot{:,kk} = t{:,kk};
    alpha_plot{:,kk} = alpha{:,kk};
    cauchyCL_plot{:,kk} = cauchyCL{:,kk};
    CL_plot{:,kk} = CL{:,kk};
    CMz_plot{:,kk} = CMz{:,kk};
%     correction{:,kk} = correction{:,kk};            
    end



%% compute cycle area    

for ll = 1: floor(length(alpha{:,kk}(time_iter_offset:end))/T_period(kk))
   
    alpha_current_cycle = alpha{:,kk}((ll-1)*T_period(kk)+1+time_iter_offset:(ll-1)*T_period(kk)+T_period(kk)+time_iter_offset);
    CL_current_cycle = CL{:,kk}((ll-1)*T_period(kk)+1:(ll-1)*T_period(kk)+T_period(kk));

    [~,idx_mincfd]  = min(alpha_current_cycle);
    [~,idx_maxcfd]  = max(alpha_current_cycle);
    
    if idx_mincfd<idx_maxcfd
        area_up(kk,ll) = trapz(alpha_current_cycle(idx_mincfd:idx_maxcfd),CL_current_cycle(idx_mincfd:idx_maxcfd));
        area_down(kk,ll) = trapz(alpha_current_cycle(idx_maxcfd:end),CL_current_cycle(idx_maxcfd:end))+trapz(alpha_current_cycle(1:idx_mincfd),CL_current_cycle(1:idx_mincfd));
    
    else
        area_down(kk,ll) = trapz(alpha_current_cycle(idx_maxcfd:idx_mincfd),CL_current_cycle(idx_maxcfd:idx_mincfd));
        area_up(kk,ll) = trapz(alpha_current_cycle(idx_mincfd:end),CL_current_cycle(idx_mincfd:end))+trapz(alpha_current_cycle(1:idx_maxcfd),CL_current_cycle(1:idx_maxcfd));
    end
end

%% compute error
% m = (CL{:,kk}(idx_maxcfd)-dynamicCL_Zanotti_CFD.CL()-(CL{:,kk}(idx_mincfd)-dynamicCL_Zanotti_CFD.CL(idx_minzano)))/(19);
% y = @(x) m*x + (CL{:,kk}(idx_mincfd)-dynamicCL_Zanotti_CFD.CL(idx_minzano));
% 
% corr = y(alpha{:,kk});
% correction{:,kk} = corr';

cd("../../")
end

%% compute error on area cycle
area = area_up-area_down;
area_err_perc_zanoEXP = (area - area_zanoEXP)/area_zanoEXP*100;
area_err_perc_zanoCFD = (area - area_zanoCFD)/area_zanoEXP*100;

table(area, area - area_zanoCFD, area - area_zanoEXP, area_err_perc_zanoCFD, area_err_perc_zanoEXP,'VariableNames',{'area','w.r.t. CFD','w.r.t. EXP','percent CFD','percent EXP'})
area_percentual_increment = (area(:,2)-area(:,1))./area(:,1)*100

%% plot CL vs alpha of the dynamic simulation

CL_alpha = figure('Position',[100,50,800,600]);
hold on;
% plot(alpha_plot{:,1},CL_plot{:,1},'color','#FFA500',"DisplayName",replace(mainFolder(1),'_',' '))
% plot(alpha_plot{:,2},CL_plot{:,2},'color','#C100FF',"DisplayName",replace(mainFolder(2),'_',' '))
plot(alpha_plot{:,3},CL_plot{:,3},'color','r',"DisplayName",replace(mainFolder(3),'_',' '))
plot(alpha_plot{:,4},CL_plot{:,4},'color','g',"DisplayName",replace(mainFolder(4),'_',' '))
% plot(alpha_plot{:,5},CL_plot{:,5},'color','c',"DisplayName",replace(mainFolder(5),'_',' '))
plot(alpha_plot{:,6},CL_plot{:,6},'color','c',"DisplayName",replace(mainFolder(6),'_',' '))

plot(dynamicCL_Zanotti_exp.alpha,dynamicCL_Zanotti_exp.CL, 'k-.',"DisplayName",'Zanotti et Al. experiment')
plot(dynamicCL_Zanotti_CFD.alpha,dynamicCL_Zanotti_CFD.CL, 'b--',"DisplayName",'Zanotti et Al. CFD')

xlabel('\alpha [°]')
ylabel('CL [-]')
title('Dynamic CL-alpha')
legend

%% plot CMz vs alpha of the dynamic simulation

CL_alpha = figure('Position',[100,50,800,600]);
hold on;
% plot(alpha_plot{:,1},CMz_plot{:,1},'color','#FFA500',"DisplayName",replace(mainFolder(1),'_',' '))
% plot(alpha_plot{:,2},CMz_plot{:,2},'color','#C100FF',"DisplayName",replace(mainFolder(2),'_',' '))
plot(alpha_plot{:,3},-CMz_plot{:,3},'color','r',"DisplayName",replace(mainFolder(3),'_',' '))
plot(alpha_plot{:,4},-CMz_plot{:,4},'color','g',"DisplayName",replace(mainFolder(4),'_',' '))
% plot(alpha_plot{:,5},CMz_plot{:,5},'color','c',"DisplayName",replace(mainFolder(5),'_',' '))
plot(alpha_plot{:,6},-CMz_plot{:,6},'color','c',"DisplayName",replace(mainFolder(6),'_',' '))

plot(dynamicCM_Zanotti_exp.alpha,dynamicCM_Zanotti_exp.CM, 'k-.',"DisplayName",'Zanotti et Al. experiment')
plot(dynamicCM_Zanotti_CFD.alpha,dynamicCM_Zanotti_CFD.CM, 'b--',"DisplayName",'Zanotti et Al. CFD')


xlabel('\alpha [°]')
ylabel('CMz [-]')
title('Dynamic CMz-alpha')
legend


%% plot integral over cycles:

figure('Position',[100,50,1000,600])
for i = 1:size(area,1)
    subplot(2,size(area,1)/2,i)
    plot(area(i,:),'o-','DisplayName','SU2')
    hold on;
    yline(area_zanoEXP,'r--','DisplayName','Zanotti exp')
    title(replace(mainFolder(i),'_',' '))
    legend
end



return
%% plot convergence criteria

figure('Position',[200,180,800,600])
subplot(1,2,1)
hold on;
% plot(time_iter{:,1}*time_step(1), cauchyCL{:,1},'color','#FFA500',"DisplayName",replace(mainFolder(1),'_',' '))
% plot(time_iter{:,2}*time_step(2), cauchyCL{:,2},'color','#C100FF',"DisplayName",replace(mainFolder(2),'_',' '))
plot(time_iter{:,3}*time_step(3), cauchyCL{:,3},'color','r',"DisplayName",replace(mainFolder(3),'_',' '))
plot(time_iter{:,4}*time_step(4), cauchyCL{:,4},'color','g',"DisplayName",replace(mainFolder(4),'_',' '))
plot(time_iter{:,5}*time_step(5), cauchyCL{:,5},'color','c',"DisplayName",replace(mainFolder(5),'_',' '))

xlabel('Simulation time')
ylabel('cauhcy CL')
title('cauchy CL over time')
legend

subplot(1,2,2)
hold on;
% plot(time_iter{:,1}*time_step(1), inner_iter{:,1},'color','#FFA500',"DisplayName",replace(mainFolder(1),'_',' '))
% plot(time_iter{:,2}*time_step(2), inner_iter{:,2},'color','#C100FF',"DisplayName",replace(mainFolder(2),'_',' '))
plot(time_iter{:,3}*time_step(3), inner_iter{:,3},'color','r',"DisplayName",replace(mainFolder(3),'_',' '))
plot(time_iter{:,4}*time_step(4), inner_iter{:,4},'color','g',"DisplayName",replace(mainFolder(4),'_',' '))
plot(time_iter{:,5}*time_step(5), inner_iter{:,5},'color','c',"DisplayName",replace(mainFolder(5),'_',' '))

xlabel('Simulation time')
ylabel('Inner iterations')
title('Inner iterations over time')
legend

%%
figure('Position',[200,180,800,500])
subplot(1,2,1)
yyaxis left
plot(time_iter{:,3}*time_step(3), cauchyCL{:,3},'color','r',"DisplayName","cauchy"+replace(mainFolder(3),'_',' '))
yyaxis right
plot(time_iter{:,3}*time_step(3), inner_iter{:,3},'color','g',"DisplayName","iter"+replace(mainFolder(3),'_',' '))
legend
title(replace(mainFolder(3),'_',' '))

subplot(1,2,2)
yyaxis left
plot(time_iter{:,4}*time_step(4), cauchyCL{:,4},'color','r',"DisplayName","cauchy"+replace(mainFolder(4),'_',' '))
yyaxis right
plot(time_iter{:,4}*time_step(4), inner_iter{:,4},'color','g',"DisplayName","iter"+replace(mainFolder(4),'_',' '))
legend
title(replace(mainFolder(4),'_',' '))


return

%% save for report .mat
cd(simulationsFolderPath+"../../Report/REPORT_RESULTS")


DYN_alpha_TC4_LM_FIRST = alpha_plot{:,6};
DYN_CL_TC4_LM_FIRST = CL_plot{:,6};
DYN_CMz_TC4_LM_FIRST = CMz_plot{:,6};


save("DYN_FIRST", ...
    "DYN_CMz_TC4_LM_FIRST","DYN_CL_TC4_LM_FIRST","DYN_alpha_TC4_LM_FIRST");