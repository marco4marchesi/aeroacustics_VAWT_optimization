function flow = acousticPostProcess(folder, dt)

%{
HELP:
function that retrieves all the acoustic interesting properties and saves
them in the struct flow
%}


% define some useful data
dt_ref = 8e-4; % [s] time for a rotation of one degree (circa) which we used as baseline value
fs = 1/dt; % sample frequency
L = 1.5; % [m] diameter
U = 8; % [m/s] , asymptotic flow velocity

N_iter_per_round = dt_ref/dt * 368;
% define transformation to decibel
spl = @(x) 20*log10(x/(20e-6));

% define rotation speed
omega=21.33;
sound_speed = sqrt(1.4*287*288.15);

% define quantities for the psd
resolutonReqd = 0.1; % [Hz]
NFFT = fs / resolutonReqd;
NFFT = 2^nextpow2(NFFT);
overlap = 0.5;

% loop on different radius length
obs_folders = dir(folder+"\observer_*");

% initialise vectors
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
    angles = [1,31];
    for k = 1:length(angles)
        ppName = "pp"+compose("%03d",angles(k))+"r"+num2str(radiuses{i});
        flow.(ppName) = readmatrix(folder+"/"+obs_folders(i).name+"/pp_FWH_"+compose("%03d",angles(k))+"_Zone_0");
        last_lap_indexes = [length(flow.(ppName))-N_iter_per_round+1,length(flow.(ppName))];
        
        % fft computation
        fftName = ppName+"_fft";
        fftName_f = ppName+"_fft_f";
        fftName_St = ppName + "_fft_St";
        fftThickness = ppName + "_fftThick";
        fftLoading = ppName + "_fftLoad";
        [flow.(fftName_f),flow.(fftName)] = fourierSingleSided(fs, flow.(ppName)(last_lap_indexes(1):last_lap_indexes(2),3)-mean(flow.(ppName)(last_lap_indexes(1):last_lap_indexes(2),3)));
        [~,flow.(fftThickness)] = fourierSingleSided(fs, flow.(ppName)(last_lap_indexes(1):last_lap_indexes(2),4)-mean(flow.(ppName)(last_lap_indexes(1):last_lap_indexes(2),4)));
        [~,flow.(fftLoading)] = fourierSingleSided(fs, flow.(ppName)(last_lap_indexes(1):last_lap_indexes(2),5)-mean(flow.(ppName)(last_lap_indexes(1):last_lap_indexes(2),5)));
        flow.(fftName_St)= flow.(fftName_f).*(L/U);
        
        % PSD computation
        psdName = ppName + "_PSD";
        psdName_f = ppName + "_PSD_f";
        psdName_St = ppName + "_PSD_St";
        window = length(flow.(ppName))/2;
        [flow.(psdName),flow.(psdName_f)] = pwelch(flow.(ppName)(last_lap_indexes(1):last_lap_indexes(2),3)-mean(flow.(ppName)(last_lap_indexes(1):last_lap_indexes(2),3)), window, overlap, NFFT, fs);
        flow.(psdName_St)= flow.(psdName_f).*(L/U);
    end
       
    
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

if history(10:100,1)~=[9:99]' % this condition happens only for unsteady simulations


flow.CD = history(end-(N_iter_per_round-1):end,9);
flow.CL = history(end-(N_iter_per_round-1):end,10);
flow.CMz = history(end-(N_iter_per_round-1):end,11);

flow.CL_phase = wrapTo2Pi(omega*dt*length(history(:,1)))+[1:N_iter_per_round]*omega*dt+pi/2;
flow.time = dt*length(history(:,1)) + ([1:N_iter_per_round]-N_iter_per_round)*dt;

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

flow.CF = sqrt(flow.CL.^2 + flow.CD.^2);
flow.dCF = diff(flow.CF)/dt;

% fft CL
[flow.CL_f,flow.CL_fft] = fourierSingleSided(fs, flow.CL-mean(flow.CL));
[flow.CF_f,flow.CF_fft] = fourierSingleSided(fs, flow.CF-mean(flow.CF));
%%
chosen_obs = [120 10]; % 1° numero è quello dell'osservatore, 2° è il raggio


X_obs = [chosen_obs(2)*cos((chosen_obs(1)-1)*2*pi/120) chosen_obs(2)*sin((chosen_obs(1)-1)*2*pi/120)];

X_prof = [-0.75*sin(flow.CL_phase') 0.75*cos(flow.CL_phase')];

time_delay = vecnorm(X_prof - repmat(X_obs,size(X_prof,1),1),2,2)/sound_speed;

flow.retarded_time = flow.time + time_delay';

flow.equivalent_speed = 0.024*sound_speed - 21.33*0.75*cos(21.33*flow.time);

end
