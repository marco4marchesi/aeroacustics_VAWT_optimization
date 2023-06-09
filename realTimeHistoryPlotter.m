%{

Matlab code for extracting data from history files real time while running
simulations

--------------------------------------------------------------------------
Author: Marco Marchesi
--------------------------------------------------------------------------

+   this code aims to extract the CD, CL, and CM values from the simulations
    and to automatically save them in a vector that contains the number of
    elements in the mesh

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

+   set

%}





%% select user

% user: set who is running the code so that the folder is chosen:
user = "doppio fisso"; % choices: "doppio fisso" "luca" ...

switch user 
    case "doppio fisso"
    matlabCodesPath = "C:\Users\marco\Desktop\UNI\2 MAGISTRALE\AEROACOUSTICS\aeroacustics_VAWT_optimization";
    %     simulationsFolderPath = "C:\Users\marco\OneDrive - Politecnico di Milano\MAGISTRALE\TerzoSemestre\CFD\PROGETTO CFD DRIVE\SIMULATIONS DRIVE\Simulations"; % drive
    simulationsFolderPath = "\\wsl.localhost\Ubuntu-20.04\home\marco\RRF_270";
end




%% init
clearvars -except matlabCodesPath simulationsFolderPath;
clc; close all;

warning('off');
% add matlab functions to the path
% rmpath(matlabCodesPath+"/polar_plotter")
% rmpath(matlabCodesPath+"/farfield_analysis")
% rmpath(matlabCodesPath+"/convergence_analysis")
% rmpath(matlabCodesPath+"/dynamic_lift_coefficient")
% addpath(matlabCodesPath)
% addpath(matlabCodesPath+"/utilitiesFunctions")

% move to simulation folder
cd(simulationsFolderPath)
% matlab_graphics;






%% ------------------------------------ CHOSE SIMULATION (FOLDER) -------------------------------------- %%
mainFolder =''; %'caseG4_ROE_highMach/';
simuFolder =''; %'caseG4_A10_newTVD\cfdG4'; % use single apices because otherwise the erase function does not work as I want
fileName = "history.dat";%'history_G4_O2_10g_newTVD.csv';%

%% LOGARITMIC PLOT

% testcase = strcat(mainFolder,simuFolder);
% cd(testcase)
R = 0.75;
c = 0.075;

try 
%% ---------------------------------------- select field ----------------------------------------------- %%
fields = ["rms_Rho_","Cauchy_CMz_","MaxCFL","CMz"]; % please, rms always first and cauchy always second!
target = [1e-8;1e-8];
zoomIter = 2000; % zoom the CD and CL plots with reference to max and min values  over the last zoomIter iterations
n_iter = 0;%80000;


stat= true;
retry = 0;
figure('Name',"Logaritmic convergence real time",'Position',[100,50,1000,600],'DockControls','on')
while(stat==true)


    %% -------------------------------------- Loop on folder simulations! Let's goooooooooooooooooooo -------------------------%%
    currentHistory = readtable(fileName);
 
    for idx_field = 1: length(fields)


        %% define from which files do the extraction


        evolution.(fields(idx_field)) = [currentHistory.(fields(idx_field))];
        n_iter_step = length(evolution.(fields(1)))-n_iter;

        if contains(fields(idx_field),"rms")
            evolution.(fields(idx_field)) = 10.^evolution.(fields(idx_field));
        end

        if contains(fields(idx_field),"CMz")
            evolution.(fields(idx_field)) = evolution.(fields(idx_field))/(R*c);
        end

        %% plot
        subplot(2,2,idx_field)
        try
            % tries to plot the updates from previous plot
            if fields(idx_field) == "CD" || fields(idx_field) == "CL" || fields(idx_field) == "CMz"
                plot(n_iter+1:n_iter+n_iter_step,evolution.(fields(idx_field))(n_iter+1:n_iter+n_iter_step),'k-','DisplayName',replace(simuFolder,'_',' '));
%                 ylim([min(evolution.(fields(idx_field))(end-min(length(evolution.(fields(idx_field)))-1,zoomIter):end))*0.99,max(evolution.(fields(idx_field))(end-min(length(evolution.(fields(idx_field)))-1,zoomIter):end))*1.01])
            else
                semilogy(n_iter+1:n_iter+n_iter_step,evolution.(fields(idx_field))(n_iter+1:n_iter+n_iter_step),'k-','DisplayName',replace(simuFolder,'_',' '));
                yline(target(idx_field),'r--','DisplayName','Convergence target')
            end


            retry = 0;
        catch

            % check if the simulation has ended. Retry one time, if good: ok. if not: tell to stop
            if retry < 5
                retry = retry+1;
            else
                stat = false;
            end
        end
        hold on;
        xlabel('Iterations')
        ylabel(replace(fields(idx_field),"_"," "))

    title(replace(fields(idx_field),"_",""))
    end

    n_iter = n_iter+n_iter_step;
    drawnow
    sgtitle("Real-time logarithmic plot of "+ simuFolder)

    %% return the result of the simulation (if ended)
    if stat==false
        fprintf('The simulation has ended! \n\n')
        if evolution.(fields(1))(end) < target(1) && evolution.(fields(2))(end) <target(2)
            fprintf('Status: convergence reached \n')
        else
            fprintf('Status: diverged \n')
        end
    else
        % wait updates from the simulation
  
        pause(2)
    end




end
catch
    cd(matlabCodesPath)

end