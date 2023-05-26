function [p_rms, thickness_rms, loading_rms, Term_1_rms, Term_2_rms, Term_3_rms, Term_4_rms]  = directivity_plot(folder,n_obs,N_laps)

Term_1_rms = [];
Term_2_rms = [];
Term_3_rms = [];
Term_4_rms = [];
thickness_rms = [];
loading_rms = [];
p_rms = [];

for i=[1:n_obs]
    
    FWH = readtable(folder+"/pp_FWH_"+compose("%03d",i)+"_Zone_0");
    
    iter_i = 1;%size(FWH,1)-4*368*N_laps;

    thickness_rms = [thickness_rms, rms(FWH.Thickness(iter_i:end))];

    loading_rms = [loading_rms, rms(FWH.Loading(iter_i:end))];

    Term_1_rms = [Term_1_rms, rms(FWH.Term_1(iter_i:end))];

    Term_2_rms = [Term_2_rms, rms(FWH.Term_2(iter_i:end))];

    Term_3_rms = [Term_3_rms, rms(FWH.Term_3(iter_i:end))];

    Term_4_rms = [Term_4_rms, rms(FWH.Term_4(iter_i:end))];

    p_rms = [p_rms, rms(FWH.P_Fluctuation(iter_i:end))];

end