function flow = acousticPostProcess(folder, dt, N_iter_per_round)

%{
HELP:
function that retrieves all the acoustic interesting properties and saves
them in the struct flow
%}

% define transformation to decibel
spl = @(x) 20*log10(x/(20e-6));

% define rotation speed
omega=21.33;

% loop on different radius length
obs_folders = dir(folder+"\observer_*");
flow.pointMatX = [];
flow.pointMatY = [];
flow.p_rmsMat = [];

for i = 1:length(obs_folders)
    n_obs(i) = length(dir(folder+"/"+obs_folders(i).name+"/pp_FWH_*"));
    
    
    % retrieve acoustic quantities
    [flow.p_rms{i}, flow.Thickness_rms{i}, flow.Loading_rms{i}, flow.Term_1{i}, flow.Term_2{i}, flow.Term_3{i}, flow.Term_4{i}]  = directivity_plot(folder+"/"+obs_folders(i).name,n_obs(i),1);
    
    % retrieve positions and angles
    radiuses{i} = str2double(erase(obs_folders(i).name,"observer_"));
    flow.points{i} = readmatrix(folder+"/"+obs_folders(i).name+"/Observer_Locations.dat");
    flow.pointMatX = [flow.pointMatX, [flow.points{i}(:,1);flow.points{i}(1,1)]];
    flow.pointMatY = [flow.pointMatY, [flow.points{i}(:,2);flow.points{i}(1,2)]];
    

    flow.theta_obs{i} = atan2(flow.points{i}(:,2),flow.points{i}(:,1));
    
    % retrieve decibel values of some quantities
    flow.SPL_p{i} = spl(flow.p_rms{i});
    flow.SPL_thickness{i} = spl(flow.Thickness_rms{i});
    flow.SPL_loading{i} = spl(flow.Loading_rms{i});
    
    flow.p_rmsMat = [flow.p_rmsMat, [flow.SPL_p{i},flow.SPL_p{i}(1)]'];

    % retrieve average noise
    flow.dB_total{i} = load(folder+"/"+obs_folders(i).name+"/pp_FWH");
    
    % recall fwh in time domain
    for k = 1:n_obs(i)
        ppName = "pp"+compose("%03d",k)+"r"+num2str(radiuses{i});
        flow.(ppName){i} = readmatrix(folder+"/"+obs_folders(i).name+"/pp_FWH_"+compose("%03d",k)+"_Zone_0");
    end
    
    last_lap_indexes = length(flow.(ppName){i})-4*360:length(flow.(ppName){i});
    ppName_f = ppName+"_f";
    ppName_fft = ppName+"_fft";
    [flow.(ppName_f){i},flow.(ppName_fft){i}] = fourierSingleSided(1/(2e-4), flow.(ppName){i}(last_lap_indexes,3)-mean(flow.(ppName){i}(last_lap_indexes,3)));
end


% retrieve aerodynamic coefficients over cycle


% ----------------------------warning-----------------------%
% mettere che controlla se ci sono history restart e prende l'ultimo, se no
% è un casotto. e poi come al solito deve prendere solo un giro se no non
% sono comparabili
% history_files = dir(folder+"\history*");
% history = readmatrix(folder+"/"+history_files(end).name);


history_files = dir(folder+"\history*");
    len = 0;
    history = [];
    if length(history_files)~=1
        for i = 1:length(history_files)
    
            h_p = readmatrix(folder+"/"+history_files(i).name);
            if folder+"/"+history_files(i).name == "E:\UNI - fisso\aeroacustica\unsteady_3_profili_deform6\history_02334.dat"
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

                iterations_lost = str2double(erase(history_files(end).name,["history_",".dat"]))+length(h_p)-len;
                if iterations_lost  > 0
                    time_loss = (iterations_lost-368*8e-4/dt) * dt ;  
                else
                    time_loss = 0;
                end
            end
        end
    else
           
        history = readmatrix(folder+"/"+history_files.name);
        iterations_lost = 0;
        time_loss = 0;
    end







flow.CD = history(end-(N_iter_per_round-1):end,9);
flow.CL = history(end-(N_iter_per_round-1):end,10);
flow.CMz = history(end-(N_iter_per_round-1):end,11);

flow.CL_phase = wrapTo2Pi(omega*dt*length(history(:,1)))+[1:N_iter_per_round]*omega*dt+pi/2;

% retrieve derivatives of the aerodynamic coefficients over cycle
flow.dCD = diff(flow.CD)/dt;
flow.dCL = diff(flow.CL)/dt;
flow.dCMz = diff(flow.CMz)/dt;
flow.dCL_phase = wrapTo2Pi(omega*dt*length(history(:,1)))+[2:N_iter_per_round]*omega*dt+pi/2;


% rms varie
flow.CD_RMS = rms(flow.CD);
flow.dCD_RMS = rms(flow.dCD);
flow.CL_RMS = rms(flow.CL);
flow.dCL_RMS = rms(flow.dCL);
flow.CMz_RMS = rms(flow.CMz);
flow.dCMz_RMS = rms(flow.dCMz);






%%format long
%%error=[100.*(flow.SPL-flowD.SPL)./(flow.SPL)];
%%media=mean(abs(error));
%%
% obs1=readtable(folder+"/pp_FWH_006_Zone_0");
% FS=10;
% figure
% tiledlayout(1,2)
% nexttile
% hold on
% plot(obs1.Time,obs1.P_Fluctuation);
% legend("Mic 135°");
% title("Pressure Fluctuation");
% xlabel('time(s)')
% ylabel("p'")
% ax = gca;
% ax.FontSize = FS;
% 
% fs = 1/(obs1.Time(2)-obs1.Time(1));
% 
% resolutonReqd = 0.1; % [Hz]
% NFFT = fs / resolutonReqd;
% NFFT = 2^nextpow2(NFFT);
% window = length(obs1.Time)/2;
% overlap = 0;
% 
% 
% [pxx1, f1] = pwelch(obs1.P_Fluctuation, window,overlap,NFFT, fs);
% nexttile
% loglog(f1, pxx1)
% legend("Mic 135°");
% title("PSD");
% xlabel('Frequency')
% ax = gca;
% ax.FontSize = FS;