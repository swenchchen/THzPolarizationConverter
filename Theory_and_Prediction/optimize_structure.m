%% 
clear all

n1=1.4:0.01:3.5;
theta_i = 15:0.1:60;
n_si = 3.418-1i*1e-7;
n_air = 1-1i*1e-7;
n_metal = 500-1i*500;

[n1_mat, theta_i_mat] = meshgrid(n1, theta_i);

for k = 1:numel(theta_i)
    thickness=0;
    rt_p = rt_trilayer(n1,n_air,n_air,thickness,theta_i(k),1,'p');
    rt_s = rt_trilayer(n1,n_air,n_air,thickness,theta_i(k),1,'s');
    ps_phase(k,:) = 180/pi*angle(rt_p.r ./ rt_s.r);
end
ps_phase(ps_phase>175) =-180;
ps_phase(ps_phase<15) =-180;
ps_phase = ps_phase+180;

theta_c = asin(1./real(n1))*180/pi;

freq = 0.2:0.01:3.5;
h=(0:0.05:60)*1e-6;
for k = 1:numel(h)
    rt_p = rt_trilayer(n_si,n_air,n_metal,h(k),30,freq,'p');
    rt_s = rt_trilayer(n_si,n_air,n_metal,h(k),30,freq,'s');
    ps_phase_h(k,:) = 180/pi*angle(rt_p.r ./ rt_s.r);
end

colorset = ColorGradient(3, 'blue');
figure('units', 'normalized', 'position', [0.1,0.05,0.2,0.45])

h1 = subplot(2,1,1);
% set(h1, 'OuterPosition', [0.026,0.7,0.95,0.2])
f_ind = interp1(freq, 1:numel(freq), [1,2,3], 'neareast');
phase_vs_h = ps_phase_h(:,f_ind);
phase_vs_h(phase_vs_h>175) = -180;
plot(h*1e6, phase_vs_h(:,1), 'linewidth', 2,'color',colorset(1,:)); hold all
plot(h*1e6, phase_vs_h(:,2), 'linewidth', 2,'color',colorset(2,:))
plot(h*1e6, phase_vs_h(:,3), 'linewidth', 2,'color',colorset(3,:))
set(h1, 'fontsize', 20); xlim([0,40]); ylim([-220,130]);
xlabel('h (\mum)'); ylabel('\phi_{ps} (бу)')
legend({'1 THz', '2 THz', '3 THz'}, 'box', 'off','location','southeast')

h2 = subplot(2,1,2);
% set(h2, 'OuterPosition', [0.06,0.03,0.9,0.26])
V=0:0.01:10;
delta_n = exp(-V);
delta_n(V<1)=exp(-1);
delta_n = delta_n/exp(-1);
plot(V, delta_n, 'linewidth', 2,'color',colorset(2,:))
xlabel('Bias'); ylabel('\Deltan')
xticks(''); yticks([0,1]); yticklabels({'min','max'})
set(h2, 'fontsize', 20);

%% SI note 1
clear all

n1=1.4:0.01:3.5;
theta_i = 15:0.1:60;
n_si = 3.418-1i*1e-7;
n_air = 1-1i*1e-7;
n_metal = 500-1i*500;
n_top = 1.4:0.01:3.5;

[n1_mat, theta_i_mat] = meshgrid(n1, theta_i);

for k = 1:numel(theta_i)
    thickness=0;
    rt_p = rt_trilayer(n1,n_air,n_air,thickness,theta_i(k),1,'p');
    rt_s = rt_trilayer(n1,n_air,n_air,thickness,theta_i(k),1,'s');
    ps_phase(k,:) = 180/pi*angle(rt_p.r ./ rt_s.r);
    
    rt_p_poly = rt_trilayer(n1,1.5-1i*1e-2,1.5-1i*1e-2,thickness,theta_i(k),1,'p');
    rt_s_poly = rt_trilayer(n1,1.5-1i*1e-2,1.5-1i*1e-2,thickness,theta_i(k),1,'s');
    ps_poly(k,:) = 180/pi*angle(rt_p_poly.r ./ rt_s_poly.r);
    
    rt_p_qz = rt_trilayer(n1,2.12-1i*1e-4,2.12-1i*1e-4,thickness,theta_i(k),1,'p');
    rt_s_qz = rt_trilayer(n1,2.12-1i*1e-4,2.12-1i*1e-4,thickness,theta_i(k),1,'s');
    ps_qz(k,:) = 180/pi*angle(rt_p_qz.r ./ rt_s_qz.r);
end
ps_phase(ps_phase>175) =-180;
ps_phase(ps_phase<15) =-180;

ps_metal = -180;
ps_phase_diff1 = ps_phase-ps_metal;

theta_c = asin(1./real(n1))*180/pi;

colorset = ColorGradient(255, 'contrast');
figure('units', 'normalized', 'position', [0.05,0.2,0.3,0.5])
h1 = subplot(2,1,2);
set(h1, 'OuterPosition', [0.05,0.03,0.85,0.55])
imagesc(theta_i, real(n1), abs(ps_phase_diff1'));  hold all
caxis([190,300]); colormap(colorset)
c = colorbar;  c.Location = 'eastoutside';
% plot(theta_c, real(n1), '-', 'color', [0.4,0.4,0.4], 'linewidth',2);
plot(theta_i, linspace(3.418,3.418,numel(theta_i)), '--', 'color', [0,0,0], 'linewidth', 1)
set(h1, 'fontsize', 20)
xlabel('\theta_i'); ylabel('n_i')

n_ind = find(abs(n1-3.418) == min(abs(n1-3.418)));
theta_ind = find(abs(theta_i-30) == min(abs(theta_i-30)));
ps_phase_si = ps_phase_diff1(:, n_ind);
h2 = subplot(2,1,1);
set(h2, 'OuterPosition', [0.04,0.56,0.733,0.2])
plot(theta_i, ps_phase_si, 'linewidth', 2)
hold all
plot(theta_i(theta_ind), ps_phase_si(theta_ind), 'ro', 'linewidth', 1,'markerfacecolor', 'w','markersize', 8)
set(h2, 'fontsize', 20)
xlim([15,60]); xticklabels([]);
ylim([230,300]);  ylabel('\Delta\phi')