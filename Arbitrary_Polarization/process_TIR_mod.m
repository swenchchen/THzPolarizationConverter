clear all

% change the name accordingly to process data measured at different bias
% duty cycles
% The number after 'air' indicates the group number
% Group 1-1 refers to LC3
% Group 1-2 refers to LC7
% Group 2-1 refers to LC6
% Group 2-2 refers to LC2
% Group 3-1 refers to LC1
% Group 3-2 refers to LC5

ref_file = 'air1.mat'; 
sam_file = 'air1_Zscan_LC1.mat';  
sam_file_woLC = 'air1_Zscan.mat'; 
save_file = 'LC1.mat';

dist = load([sam_file_woLC(1:end-4),'_distance.txt']);
air_gap = (dist(1,:)-dist(1,1))*1e-6+5.7e-6;
target_phase = -127.5;
freq_range = [0.6,2.4];
win_mode = 2;
win_length = 3.4;
ref_load = process_terasmart(ref_file,0,0,win_mode,0,win_length);
sam_load = process_terasmart(sam_file,0,0,win_mode,0,win_length);
sam_woLC_load = process_terasmart(sam_file_woLC,0,0,win_mode,0,win_length);
freq = ref_load.data_fd(1,:);
time = ref_load.data_itp(1,:);
sam_p_fd = sam_load.data_fd(2:12,:);
sam_s_fd = sam_load.data_fd(13:23,:);
sam_woLC_p_fd = sam_woLC_load.data_fd(2:12,:);
sam_woLC_s_fd = sam_woLC_load.data_fd(13:23,:);
ref_p_fd = ref_load.data_fd(2,:);
ref_s_fd = ref_load.data_fd(3,:);
ref_ps_exp = ref_p_fd./ref_s_fd;

n_si = 3.419-1i*1e-4;
n_air = 1-1i*1e-2;
metal_sigma = 4*10^7;
d_metal = 200e-9;
metal_prop = nk_epsilon_cond(metal_sigma*d_metal*1000,freq,d_metal,1,'conductivity');
n_metal = metal_prop.nk;
theta_i=30;

rtp_ref = rt_trilayer(n_si,n_air,n_air,0,theta_i,freq,'p');
rts_ref = rt_trilayer(n_si,n_air,n_air,0,theta_i,freq,'s');
ref_ps_ratio = rtp_ref.r ./ rts_ref.r;

exp_ratio_norm = (sam_p_fd ./ sam_s_fd) ./ repmat(ref_ps_exp, 11, 1);
exp_ratio_cal = exp_ratio_norm.*repmat(ref_ps_ratio, 11, 1); 
exp_woLC_ratio_norm = (sam_woLC_p_fd ./ sam_woLC_s_fd) ./ repmat(ref_ps_exp, 11, 1);
exp_woLC_ratio_cal = exp_woLC_ratio_norm.*repmat(ref_ps_ratio, 11, 1); 
colorset = ColorGradient(11, 'full');
figure
for k = 1:11
    rtp_gap = rt_trilayer(n_si,n_air,n_metal,air_gap(k),theta_i,freq,'p');
    rts_gap = rt_trilayer(n_si,n_air,n_metal,air_gap(k),theta_i,freq,'s');
    gap_ps_ratio = rtp_gap.r ./ rts_gap.r;
    gap_phase(k,:) = unwrap(angle(gap_ps_ratio))*180/pi;
    
    subplot(1,3,1)
    plot(freq, 180/pi*unwrap(angle(exp_ratio_cal(k,:))), 'color', colorset(k,:), 'linewidth',1); hold all
    plot(freq, 180/pi*unwrap(angle(exp_woLC_ratio_cal(k,:))), '--', 'color', colorset(k,:), 'linewidth',1)
    plot(freq, gap_phase(k,:), 'k')
    subplot(1,3,2)
    plot(freq, 180/pi*unwrap(angle(exp_ratio_cal(k,:))-angle(exp_woLC_ratio_cal(k,:))), 'color', colorset(k,:), 'linewidth',1)
    hold all; xlim([0.5,3.5])
    subplot(1,3,3)
    plot(freq, abs(exp_ratio_cal(k,:)), 'color', colorset(k,:), 'linewidth',1); hold all; xlim([0.5,3.5])
    hold all
    plot(freq, abs(exp_woLC_ratio_cal(k,:)),'--', 'color', colorset(k,:), 'linewidth',1); hold all; xlim([0.5,3.5])
end
subplot(1,3,1)
plot(freq, linspace(target_phase,target_phase,numel(freq)), 'k--')
xlim(freq_range); grid('on')

save(save_file, 'freq', 'exp_ratio_cal','exp_woLC_ratio_cal');
