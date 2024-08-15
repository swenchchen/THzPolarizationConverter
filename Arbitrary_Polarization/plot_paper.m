clear all

target_phase = [85, -180, 52.5, 0, -95, 0, -127.5];
target_amp = [1.1636, 0.01, 1.9597, 0, -1.1636, -1.7321, -1.9597];
target_h = [47.02, 1.13, 30.58, 11.64, 6.06, 22.14, 5.9]*1e-6;
f_down = [0.64, 0.2, 0.71, 1.215, 1, 0.64, 0.67];
f_up = [3.365, 3.5, 2.33, 3.38, 3, 1.78, 2.34];
flip_ind = [1, 1, 1, 1, 1, 1, 1];
best_ind = [6,2,6,2,4,6,5]; % index that fits the target phase best, i.e. the best height

n_si = 3.418-1i*1e-4;
n_air = 1-1i*1e-4;
n_metal = 500-1i*500;
theta_i = 30;

figure
for k = 1:7
    file_name = ['LC',int2str(k),'.mat'];
    load(file_name)
    phase_exp = zeros(11, numel(freq));
    fit_ind = find(freq>0.5 & freq<3);
    for d = 1:11
        phase_exp(d,:) = (unwrap(angle(exp_ratio_cal(d,:))+pi/5)-pi/5)*180/pi;
    end
    compare_ind = find(freq > f_down(k) & freq<f_up(k));

    data_out(k).freq = freq(compare_ind);
    data_out(k).phase = flip_ind(k)*phase_exp(best_ind(k), compare_ind);
    data_out(k).amp = abs(exp_ratio_cal(best_ind(k), compare_ind));
    data_out(k).TIRphase = angle(exp_woLC_ratio_cal(best_ind(k), fit_ind));
    
    subplot(1,2,1)
    plot(data_out(k).freq, data_out(k).phase); hold all
    plot(data_out(k).freq, flip_ind(k)*linspace(target_phase(k), target_phase(k), numel(data_out(k).freq)), 'k--')
    
    % 拟合确定实测金属镜高度
    h_array = linspace(0.7*target_h(k), max([1.3*target_h(k), target_h(k)+1e-6]), 40);
    for j = 1:numel(h_array)
        rtp_gap = rt_trilayer(n_si,n_air,n_metal,h_array(j),theta_i,freq(fit_ind),'p');
        rts_gap = rt_trilayer(n_si,n_air,n_metal,h_array(j),theta_i,freq(fit_ind),'s');
        gap_ps_ratio = rtp_gap.r ./ rts_gap.r;
        TIR_err(j) = mean(abs(data_out(k).TIRphase - angle(gap_ps_ratio)).^2);
    end
    TIR_ind = find(TIR_err == min(TIR_err));
    h_best(k) = h_array(TIR_ind);
    rtp_gap = rt_trilayer(n_si,n_air,n_metal,h_best(k),theta_i,freq(fit_ind),'p');
    rts_gap = rt_trilayer(n_si,n_air,n_metal,h_best(k),theta_i,freq(fit_ind),'s');
    gap_ps_ratio = rtp_gap.r ./ rts_gap.r;
    subplot(1,2,2)
    plot(freq(fit_ind), data_out(k).TIRphase,'o'); hold all
    plot(freq(fit_ind), angle(gap_ps_ratio), 'linewidth', 2)
end
h_best([3,7,6,1,5,4,2])*1e6  % the actual mirror height values for the 6 groups

%%

plot_ind = [3,7,6,2,1,5]; 
target_amp = [1.1636, 1.7321, 1.9597, 0, -1.1636, -1.7321, -1.9597]; 

t = linspace(0,2*pi, 100);
Es = cos(t);

figure('units', 'normalized', 'position',[0.1,0.05,0.5,0.75])
for k = 1:numel(plot_ind)
    ind = plot_ind(k);
    amp_corr = data_out(ind).amp / mean(data_out(ind).amp) * target_amp(ind);
    colorset = ColorGradient(numel(data_out(ind).freq),'full');
    for j = 1:4:numel(data_out(ind).freq)
        Ep(k).E(j,:) = amp_corr(j).*cos(t+data_out(ind).phase(j)*pi/180);
        E_abs = sqrt(1+abs(amp_corr(j))^2);
        freq_fix = repmat(data_out(ind).freq(j), 1, numel(t));
        plot_pos = [1,4,2,5,3,6];
        subplot(2,3,plot_pos(k))
        plot(Ep(k).E(j,:)/E_abs,  Es/E_abs, 'color', [colorset(j,:)] ,'linewidth', 1); hold all
        xlim([-1,1]); ylim([-1,1]); daspect([1,1,1])
    end
    set_phase = target_phase(ind)*pi/180;
    set_amp = target_amp(ind);
    E_abs = sqrt(1+set_amp^2);
    Ep_theory = set_amp.*cos(t+set_phase);
    if k==3 || k==4
    else
        plot(Ep_theory/E_abs,  Es/E_abs, 'k:' ,'linewidth', 1.5)
    end
    set(gca, 'fontsize', 18); xlabel('E_p'); ylabel('E_s')
    ch = colorbar;  colormap(colorset);
    freq_range = round([data_out(ind).freq(1), data_out(ind).freq(end)]*10)/10;
    f_num = numel(data_out(ind).freq);
    freq_ticks = round([data_out(ind).freq(1), data_out(ind).freq(round(f_num/2)), data_out(ind).freq(end)]*10)/10;
    ch.Ticks = freq_ticks;
    ch.Location = 'northoutside';
    ch.FontSize = 18;
    caxis(freq_range);
end

%%
colorset2 = ColorGradient(3, '7lines');
plot_ind = [3,7,6,2,1,5];
figure('units', 'normalized', 'position',[0.1,0.1,0.32,0.65])
gray_color = [0.6,0.6,0.6];

f = [1,2,3,4];
v1 = [0.2, 42.5; 3.5, 42.5; 3.5, 62.5; 0.2, 62.5];
patch('Faces', f, 'Vertices', v1, 'FaceColor', [0.8,0.8,0.8], 'FaceAlpha',0.5,'EdgeColor', 'none'); hold all; box('on')
v2 = [0.2, -137.5; 3.5, -137.5; 3.5, -117.5; 0.2, -117.5];
patch('Faces', f, 'Vertices', v2, 'FaceColor', [0.8,0.8,0.8], 'FaceAlpha',0.5,'EdgeColor', 'none');
v3 = [0.2, -10; 3.5, -10; 3.5, 10; 0.2, 10];
patch('Faces', f, 'Vertices', v3, 'FaceColor', [0.8,0.8,0.8], 'FaceAlpha',0.5,'EdgeColor', 'none');
v4 = [0.2, -190; 3.5, -190; 3.5, -170; 0.2, -170];
patch('Faces', f, 'Vertices', v4, 'FaceColor', [0.8,0.8,0.8], 'FaceAlpha',0.5,'EdgeColor', 'none');
v5 = [0.2, 75; 3.5, 75; 3.5, 95; 0.2, 95];
patch('Faces', f, 'Vertices', v5, 'FaceColor', [0.8,0.8,0.8], 'FaceAlpha',0.5,'EdgeColor', 'none');
v6 = [0.2, -105; 3.5, -105; 3.5, -85; 0.2, -85];
patch('Faces', f, 'Vertices', v6, 'FaceColor', [0.8,0.8,0.8], 'FaceAlpha',0.5,'EdgeColor', 'none');

for k = 1:numel(plot_ind)
    ind = plot_ind(k);
    if k == 1
        load('phase_20_n20.mat');
        set_color1 = 0.9*colorset2(k,:);
        set_color2 = 0.3*(1-colorset2(k,:));
        h1=plot(phase_predict(1).data(1,:), phase_predict(1).data(2,:),'linewidth', 4, 'color', set_color1); hold all
        h2=scatter(data_out(ind).freq(1:4:end), data_out(ind).phase(1:4:end), 25,...
            'o', 'linewidth', 0.5, 'markeredgecolor', set_color2, 'markerfacecolor', 'w','markerfacealpha',0.8);
    elseif k==2
        plot(phase_predict(2).data(1,:), phase_predict(2).data(2,:),'linewidth', 4, 'color', set_color1)
        scatter(data_out(ind).freq(1:4:end), data_out(ind).phase(1:4:end), 25,...
            'o', 'linewidth', 0.5,  'markeredgecolor', set_color2, 'markerfacecolor', 'w','markerfacealpha',0.8)
    elseif k==3 
        set_color1 = colorset2(k-1,:);
        set_color2 = [0,0,0];
        load('phase_n30_0.mat');
        h3=plot(phase_predict(2).data(1,:), phase_predict(2).data(2,:),'linewidth', 4, 'color', set_color1);
        h4=scatter(data_out(ind).freq(1:4:end), data_out(ind).phase(1:4:end), 26,...
            'd', 'linewidth', 0.5, 'markeredgecolor', set_color2, 'markerfacecolor', 'w','markerfacealpha',0.8);
    elseif k==4
        plot(phase_predict(1).data(1,:), phase_predict(1).data(2,:),'linewidth', 4, 'color', set_color1)
        scatter(data_out(ind).freq(1:4:end), data_out(ind).phase(1:4:end), 26,...
            'd', 'linewidth', 0.5, 'markeredgecolor', set_color2, 'markerfacecolor', 'w','markerfacealpha',0.8)
    elseif k==5
        set_color1 = colorset2(k-2,:);
        set_color2 = [0,0,0];
        load('phase_15_n40.mat');
        h5=plot(phase_predict(1).data(1,:), phase_predict(1).data(2,:),'linewidth', 4, 'color', set_color1);
        h6=scatter(data_out(ind).freq(1:4:end), data_out(ind).phase(1:4:end), 30,...
            's', 'linewidth', 0.5, 'markeredgecolor', set_color2, 'markerfacecolor', 'w','markerfacealpha',0.8);
    else
        plot(phase_predict(2).data(1,:), phase_predict(2).data(2,:),'linewidth', 4, 'color', set_color1)
        scatter(data_out(ind).freq(1:4:end), data_out(ind).phase(1:4:end), 30,...
            's', 'linewidth', 0.5, 'markeredgecolor', set_color2, 'markerfacecolor', 'w','markerfacealpha',0.8)
    end
end
xlim([0.2,3.5]); ylim([-200,120]); %yticks([-180,-90,0,90]);
set(gca, 'fontsize', 22);
xlabel('\omega/2\pi (THz)'); ylabel('\phi_{ps} (°)')

%%
colorset2 = ColorGradient(7, '7lines');
plot_ind = [3,7,6,2,1,5]; % 将第二组（表里最后一行）改变Ep/Es幅值比，使之输出-30度线偏
figure('units', 'normalized', 'position',[0.1,0.1,0.27,0.5])

for k = 1:numel(plot_ind)
    ind = plot_ind(k);
    set_color1 = colorset2(k,:);
    plot(data_out(ind).freq(1:1:end), data_out(ind).amp(1:1:end),...
        '-', 'linewidth', 2, 'color', set_color1);
    hold all
end
xlim([0.2,3.5]); 
set(gca, 'fontsize', 20);
xlabel('\omega/2\pi (THz)'); ylabel('|r_p/r_s|')
legend({'1-1 \Deltan=0.111', '1-2 \Deltan=0.289', '2-1 \Deltan=0.241','2-2 \Deltan=0.106',...
    '3-1 \Deltan=0.037', '3-2 \Deltan=0.200'},'fontsize',16,'box', 'off')