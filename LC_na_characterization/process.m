%% 
% Calculate LC nk by single transmission & multiple transmission
% The LC thickness can be determined by the FP period
% Correct the pulse shift (induced by Si window thickness error) by matching the FP oscillations

clear all

data90_win2 = process_terasmart('LC1_90deg_0V.mat',0,0,2,0,4);
si_ref_win2 = process_terasmart('ref-Si.mat',0,0,2,0,4);

freq = data90_win2.data_fd(1,:);
sol_ind = find(freq>0.2 & freq<3.5);
freq_sol = freq(sol_ind);

ratio1 = data90_win2.data_fd(2,:) ./ si_ref_win2.data_fd(2,:);
ratio2 = data90_win2.data_fd(3,:) ./ si_ref_win2.data_fd(3,:);

grad1 = -0.1;
grad2 = 0.06;

ratio1_sol = ratio1(:, sol_ind).*exp(1i*freq_sol*grad1);
ratio2_sol = ratio2(:, sol_ind).*exp(1i*freq_sol*grad2);

n_si = 3.418-1i*1e-4;
n_air = 1-1i*1e-4;
n_LC = 1.6-1i*1e-2;
N_sam(1).n = n_si; N_sam(2).n = n_LC; N_sam(3).n = n_si; 
N_ref(1).n = n_si; N_ref(2).n = n_air; N_ref(3).n = n_si; 
solve_num = [2, 0];
d = 210e-6; d_sam = d;  d_ref = d;
theta_i = 0;
rt2(1).rt = 't_single';  rt2(2).rt = 't_single';  
ps2(1).ps = 'p';  ps2(2).ps = 'p';  
accuracy = 1e-4;

nk_result1 = nk_solver(freq_sol, ratio1_sol, N_sam, N_ref, solve_num, d_sam, d_ref, theta_i, rt2, ps2, accuracy);
nk_result2 = nk_solver(freq_sol, ratio2_sol, N_sam, N_ref, solve_num, d_sam, d_ref, theta_i, rt2, ps2, accuracy);

nk1 = nk_result1.n - 1i*nk_result1.k;
nk2 = nk_result2.n - 1i*nk_result2.k;

data90_win3 = process_terasmart('LC1_90deg_0V.mat',0,0,3,0,10);
si_ref_win3 = process_terasmart('ref-Si.mat',0,0,3,0,10);

ratio_win3_1 = data90_win3.data_fd(2,:) ./ si_ref_win3.data_fd(2,:);
ratio_win3_2 = data90_win3.data_fd(3,:) ./ si_ref_win3.data_fd(3,:);

rt_ch1 = rt_trilayer(n_si,nk1,n_si,d,theta_i,freq_sol,'p');
rt_ch2 = rt_trilayer(n_si,nk2,n_si,d,theta_i,freq_sol,'p');
rt_si = rt_trilayer(n_si,n_air,n_si,d,theta_i,freq_sol,'p');
t_ch1 = rt_ch1.t ./ rt_si.t_single;
t_ch2 = rt_ch2.t ./ rt_si.t_single;

figure
subplot(2,2,1)
plot(freq, abs(ratio_win3_1)); hold all
plot(freq_sol, abs(t_ch1)); xlim([0.2,3.5])
subplot(2,2,2)
plot(freq, abs(ratio_win3_2)); hold all
plot(freq_sol, abs(t_ch2)); xlim([0.2,3.5])
subplot(2,2,3)
plot(freq_sol, nk_result1.n); hold all
plot(freq_sol, nk_result2.n)
subplot(2,2,4)
plot(freq_sol, nk_result1.a); hold all
plot(freq_sol, nk_result2.a)

%%
save('LC_nk.mat', 'freq_sol', 'nk_result1', 'nk_result2');
colorset = ColorGradient(7, '7lines');
figure
subplot(1,2,1)
yyaxis left
plot(freq_sol, nk_result1.n, 'linewidth', 1.5,'color', colorset(1,:)); hold all
plot(freq_sol, nk_result2.n, '--', 'linewidth', 1.5,'color', colorset(1,:))
xlabel('\omega/2\pi (THz)'); ylabel('n'); ylim([1.45, 1.85])
set(gca, 'fontsize', 20)
yyaxis right
plot(freq_sol, nk_result1.n-nk_result2.n, 'linewidth', 1.5)
legend({'n_e', 'n_o', '\Deltan'}, 'edgecolor', 'none')
ylim([0.2,0.4]); ylabel('\Deltan')
subplot(1,2,2)
plot(freq_sol, nk_result1.a, 'linewidth', 1.5,'color', colorset(1,:)); hold all
plot(freq_sol, nk_result2.a, '--', 'linewidth', 1.5,'color', colorset(1,:))
xlabel('\omega/2\pi (THz)'); ylabel('\alpha (cm^{-1})')
set(gca, 'fontsize', 20)
legend({'\alpha_e', '\alpha_o'}, 'edgecolor', 'none', 'location', 'northwest')