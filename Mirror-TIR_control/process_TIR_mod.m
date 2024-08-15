%% 
clear all

ref_file = 'air1.mat'; % Empty prism without LC
sam_file = 'air1_Z_scan_full.mat';  % mirror-coupled prism without LC

win_mode = 2;
win_length = 8;
ref_load = process_terasmart(ref_file,0,0,win_mode,0,win_length);
sam_load = process_terasmart(sam_file,0,0,win_mode,0,win_length);
freq = ref_load.data_fd(1,:);
time = ref_load.data_itp(1,:);
sam_p_td = sam_load.data_itp(2:12,:);
sam_s_td = sam_load.data_itp(13:23,:);
sam_p_fd = sam_load.data_fd(2:12,:);
sam_s_fd = sam_load.data_fd(13:23,:);
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
for k=1:11
    phase_exp(k,:) = 180/pi*unwrap(angle(exp_ratio_cal(k,:)))-360;
end

dist = load('air1_Z_scan_full_distance.txt');
air_gap = (dist(1,:)-dist(1,1))*1e-6+1.4e-6;

grad_ind = find(freq>1 & freq<3);
freq_range = freq(grad_ind(end))-freq(grad_ind(1));
for k = 1:numel(air_gap)
    d_var = linspace(max([air_gap(k)-1e-6,0]), air_gap(k)+1.5e-6,50);
    for j = 1:numel(d_var)
        rtp_gap = rt_trilayer(n_si,n_air,n_metal,d_var(j),theta_i,freq,'p');
        rts_gap = rt_trilayer(n_si,n_air,n_metal,d_var(j),theta_i,freq,'s');
        gap_ps_ratio(j,:) = rtp_gap.r ./ rts_gap.r;
        phase_err(j) = mean(abs(angle(gap_ps_ratio(j,grad_ind))-angle(exp_ratio_cal(k,grad_ind))));
    end
    min_err_ind(k) = find(phase_err == min(phase_err));
    d_fit(k) = d_var(min_err_ind(k));
    gap_ps_ratio_best(k,:) = gap_ps_ratio(min_err_ind(k),:);
    gap_ps_phase(k,:) = 180/pi*angle(gap_ps_ratio_best(k,:));
    grad(k) = gap_ps_phase(k,grad_ind(end))-gap_ps_phase(k,grad_ind(1));
end


colorset = ColorGradient(11, 'full');
figure('units', 'normalized','position',[0.1,0.1,0.7,0.8])
for k = 1:11
    subplot(2,2,1)
    plot(freq(1:4:end), phase_exp(k,1:4:end), 'o', 'color', colorset(k,:), 'linewidth',1); hold all
    plot(freq, 180/pi*angle(gap_ps_ratio_best(k,:)), 'color', colorset(k,:), 'linewidth',2)
    xlim([0.3,3.5]); 
    subplot(2,2,2)
    plot(freq(1:4:end), abs(exp_ratio_cal(k,(1:4:end))), 'o', 'color', colorset(k,:), 'linewidth',1); hold all
    xlim([0.3,3.5]);
end
plot(freq, linspace(1,1,numel(freq)), '--', 'color', [0.5,0.5,0.5], 'linewidth',2)
subplot(2,2,1)
xlabel('\omega/2\pi (THz)'); ylabel('\phi_{ps} (бу)')
set(gca, 'fontsize', 20)
% legend({'set 1.40\mum', 'fit 1.67\mum', 'set 1.9\mum', 'fit 2.02\mum'})
subplot(2,2,2)
xlabel('\omega/2\pi (THz)'); ylabel('|r_p/r_s|')
set(gca, 'fontsize', 20)

subplot(2,2,3)
for k = 1:11
    plot(time-241.5, sam_p_td(k,:), 'color', colorset(k,:), 'linewidth', 0.5); hold all
    plot(time-251.5, sam_s_td(k,:)-10, 'color', colorset(k,:), 'linewidth', 0.5); hold all
end
xlim([0,4])
xlabel('Time (ps)'); ylabel('Amp. (a.u.)');
set(gca, 'fontsize', 20)

subplot(2,2,4)
yyaxis left
plot(air_gap*1e6, 'b^-', 'linewidth', 1.5, 'markerfacecolor', 'w'); hold all
plot(d_fit*1e6, 'ro--', 'linewidth', 1.5, 'markerfacecolor', 'w')
xlim([0,12])
ylabel('h (\mum)'); 
yyaxis right
plot((air_gap-d_fit)*1e6, 'v-', 'linewidth', 1, 'markerfacecolor', 'w')
xlabel('No. of measurements'); ylabel('\Deltah (\mum)')
set(gca, 'fontsize', 20); grid('on')
legend({'Set height', 'Fit height', 'Height error'}, 'location', 'northwest', 'box', 'off')