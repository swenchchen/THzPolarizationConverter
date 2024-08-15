clear all

file_num = [0.999, 0.9:-0.1:0.7, 0.5:-0.1:0.3, 0.25, 0.2, 0.08, 0.06, 0.05:-0.01:0.02, 0.015:-0.005:0.005, 0.001];

data_num = numel(file_num);
colorset = ColorGradient(data_num, 'full');
figure ('units', 'normalized', 'position', [0,0.1,0.9,0.5])
for k = 1:data_num
    sam_file = [num2str(file_num(k)), '_back_LC.mat'];
    ref_file = [num2str(file_num(k)), '_back_air.mat'];

    sam_load = process_terasmart(sam_file, 0,0,2,0,4);
    ref_load = process_terasmart(ref_file, 0,0,2,0,4);
    
    freq = sam_load.data_fd(1,:);
    time = sam_load.data_itp(1,:);
    sam_ps = sam_load.data_fd(2,:) ./ sam_load.data_fd(3,:);
    ref_ps = ref_load.data_fd(2,:) ./ ref_load.data_fd(3,:);
    ps_norm(k,:) = sam_ps./ref_ps;
    ps_phase(k) = unwrap_TDS(freq, ps_norm(k,:), [0.4,1.8],3);
    
    sam_p_td(k,:) = sam_load.data_win(2,:);
    sam_s_td(k,:) = sam_load.data_win(3,:);
    
    subplot(1,3,1)
    plot(time-249.5, sam_p_td(k,:), 'color', colorset(k,:))
    hold all
    plot(time(1:end-299)-249.5, sam_s_td(k,300:end)-6, 'color', colorset(k,:))
    hold all
    xlim([0,5])
    subplot(1,3,2)
    plot(freq, unwrap_TDS(freq, ps_norm(k,:), [0.4,1.8],1)*180/pi, 'color', colorset(k,:)); hold all
    xlim([0.2,3.5])
    subplot(1,3,3)
    plot(freq, abs(ps_norm(k,:)), 'color', colorset(k,:)); hold all
    xlim([0.2,3.5])
end
subplot(1,3,1)
xlabel('Time (ps)'); ylabel('Amplitude (a.u.)')
set(gca, 'fontsize', 20)
subplot(1,3,2)
xlabel('\omega/2\pi (THz)'); ylabel('\phi_{ps} (бу)')
set(gca, 'fontsize', 20)
subplot(1,3,3)
xlabel('\omega/2\pi (THz)'); ylabel('|t_p/t_s|')
set(gca, 'fontsize', 20)

figure
plot(100*[1,file_num(2:end-1),0], ps_phase, 'o-', 'linewidth', 1)

%%
backward_dc = 100*[1,file_num(2:end-1),0];
backward_phase = ps_phase;
save('backward_bias.mat', 'backward_dc', 'backward_phase');

load('forward_bias.mat')
d=210e-6;
phase_deltan_coef = 2*pi*1e12*d/3e8;

colorset = ColorGradient(7, '7lines');
figure('units', 'normalized', 'position', [0.1,0.1,0.8,0.8])
load('LC_nk.mat');
subplot(2,2,1)
yyaxis left
plot(freq_sol, nk_result1.n, 'linewidth', 1.5,'color', colorset(1,:)); hold all
plot(freq_sol, nk_result2.n, '--', 'linewidth', 1.5,'color', colorset(1,:))
xlabel('\omega/2\pi (THz)'); ylabel('n'); ylim([1.45, 1.85])
set(gca, 'fontsize', 20)
yyaxis right
plot(freq_sol, nk_result1.n-nk_result2.n, 'linewidth', 1.5)
legend({'n_e', 'n_o', '\Deltan'}, 'box', 'off','NumColumns',2)
xlim([0.1,3.6]); ylim([0.2,0.4]); ylabel('\Deltan')
subplot(2,2,2)
plot(freq_sol, nk_result1.a, 'linewidth', 1.5,'color', colorset(1,:)); hold all
plot(freq_sol, nk_result2.a, '--', 'linewidth', 1.5,'color', colorset(1,:))
xlabel('\omega/2\pi (THz)'); ylabel('\alpha (cm^{-1})')
xlim([0.1,3.6]); set(gca, 'fontsize', 20)
legend({'\alpha_e', '\alpha_o'}, 'box', 'off', 'location', 'northwest')

load('LC_response.mat');
subplot(2,2,3)
plot(meas_time, -delta_n, 'linewidth', 1)
xlabel('Measurement time (s)');
ylabel('\Deltan')
set(gca, 'fontsize', 20)
grid('on'); xlim([0,320])

subplot(2,2,4)
plot(forward_dc, -forward_phase/phase_deltan_coef, 'o-', 'linewidth', 1, 'markerfacecolor', 'w', 'markersize', 6)
hold all
plot(backward_dc, -backward_phase/phase_deltan_coef, '^-', 'linewidth', 1, 'markerfacecolor', 'w', 'markersize', 6)
xlabel('Duty Cycle (%)'); ylabel('\Deltan')
set(gca, 'fontsize', 20)
legend({'Forward', 'Backward'}, 'fontsize', 18,'box', 'off','location', 'southwest')
xlim([0,10]); ylim([0.05,0.35]); grid('on')

%% inset figures
figure('units', 'normalized', 'position', [0.1,0.1,0.27,0.22])
subplot(1,2,1)
plot(meas_time, -delta_n)
set(gca, 'fontsize', 20)
xlim([0,15])

subplot(1,2,2)
plot(forward_dc, -forward_phase/phase_deltan_coef, 'o-', 'linewidth', 1, 'markerfacecolor', 'w', 'markersize', 4)
hold all
plot(backward_dc, -backward_phase/phase_deltan_coef, '^-', 'linewidth', 1, 'markerfacecolor', 'w', 'markersize', 4)
set(gca, 'fontsize', 20)