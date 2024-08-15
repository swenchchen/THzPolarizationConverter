clear all

load('LC0.031.mat');
ratio_LC(1:2,:) = exp_ratio_cal([4,8],:);
ratio_woLC(1:2,:) = exp_woLC_ratio_cal([4,8],:);
load('LC0.05.mat');
ratio_LC(3,:) = exp_ratio_cal(4,:);
ratio_woLC(3,:) = exp_woLC_ratio_cal(4,:);
load('LC08.mat');
ratio_LC(4,:) = exp_ratio_cal(3,:);
ratio_woLC(4,:) = exp_woLC_ratio_cal(3,:);

%%
colorset = ColorGradient(7,'7lines');
figure('units', 'normalized', 'position', [0.1,0.1,0.25,0.7])
subplot(2,1,1)
v1 = [1, -200; 1.6, -200; 1.6, 110; 1, 110];
f = [1,2,3,4];
v2 = [3.4, -200; 3.7, -200; 3.7, 110; 3.4, 110];
patch('Faces', f, 'Vertices', v1, 'FaceColor', [0.8,0.8,0.8], 'FaceAlpha',0.5,'EdgeColor', 'none'); hold all
patch('Faces', f, 'Vertices', v2, 'FaceColor', [0.8,0.8,0.8], 'FaceAlpha',0.5,'EdgeColor', 'none')

h1=plot(freq, 180/pi*unwrap(angle(ratio_LC(1,:))), 'color', colorset(1,:), 'linewidth',2); hold all
h2=plot(freq, 180/pi*unwrap(angle(ratio_LC(2,:))), 'color', colorset(2,:), 'linewidth',2);
h3=plot(freq, 180/pi*unwrap(angle(ratio_LC(3,:))), 'color', colorset(3,:), 'linewidth',2);
h4=plot(freq, 180/pi*unwrap(angle(ratio_LC(4,:))), 'color', colorset(4,:), 'linewidth',2);

plot(freq, linspace(0,0,numel(freq)), 'k--');
plot(freq, linspace(-90,-90,numel(freq)), 'k--');
plot(freq, linspace(-180,-180,numel(freq)), 'k--');
plot(freq, linspace(90,90,numel(freq)), 'k--');
xlim([1,3.7]); ylim([-200,110]); box('on')
yticks([-180,-90 ,0, 90]);  grid('on');
xlabel('\omega/2\pi (THz)'); ylabel('\phi_{ps} (бу)')
set(gca, 'fontsize', 20)
legend([h1,h2,h3,h4],{'Pol. 1', 'Pol. 2', ...
    'Pol. 3',  'Pol. 4'},'box', 'off','NumColumns', 2) 

subplot(2,1,2)
v3 = [1, 0.7; 1.6, 0.7; 1.6, 1.2; 1, 1.2];
v4 = [3.4, 0.7; 3.7, 0.7; 3.7, 1.2; 3.4, 1.2];
patch('Faces', f, 'Vertices', v3, 'FaceColor', [0.8,0.8,0.8], 'FaceAlpha',0.5,'EdgeColor', 'none'); hold all
patch('Faces', f, 'Vertices', v4, 'FaceColor', [0.8,0.8,0.8], 'FaceAlpha',0.5,'EdgeColor', 'none')

h1 = plot(freq, abs(ratio_LC(1,:)), 'color', colorset(1,:), 'linewidth',2); hold all
h2 = plot(freq, abs(ratio_LC(2,:)), 'color', colorset(2,:), 'linewidth',2);
h3 = plot(freq, abs(ratio_LC(3,:)), 'color', colorset(3,:), 'linewidth',2);
h4 = plot(freq, abs(ratio_LC(4,:)), 'color', colorset(4,:), 'linewidth',2);
xlim([1,3.7]); ylim([0.7,1.2]); box('on'); grid('on');
xlabel('\omega/2\pi (THz)'); ylabel('|r_p/r_s|')
set(gca, 'fontsize', 20)
legend([h1,h2,h3,h4],{'Pol. 1', 'Pol. 2', ...
    'Pol. 3',  'Pol. 4'},'box', 'off','NumColumns', 2,'location','south') 
% legend([h1,h2,h3,h4],{'DC = 3.1%, h = 1.7\mum', 'DC = 3.1%, h = 4.6\mum', ...
%     'DC = 5%, h = 12.3\mum',  'DC = 80%, h = 24\mum'},'box', 'off','location','southwest') 

%%
figure('units', 'normalized', 'position', [0.1,0.1,0.4,0.5])

LC_phase(1,:) = 180/pi*unwrap_TDS(freq, ratio_LC(1,:)./ratio_woLC(1,:),[0.5,2.5],1);
LC_phase(2,:) = 180/pi*unwrap_TDS(freq, ratio_LC(2,:)./ratio_woLC(2,:),[0.5,2.5],1);
LC_phase(3,:) = 180/pi*unwrap_TDS(freq, ratio_LC(3,:)./ratio_woLC(3,:),[0.5,2.5],1);
LC_phase(4,:) = 180/pi*unwrap_TDS(freq, ratio_LC(4,:)./ratio_woLC(4,:),[0.5,2.5],1);

grad(1) = unwrap_TDS(freq, ratio_LC(1,:)./ratio_woLC(1,:),[0.5,2.5],3);
grad(2) = unwrap_TDS(freq, ratio_LC(2,:)./ratio_woLC(2,:),[0.5,2.5],3);
grad(3) = unwrap_TDS(freq, ratio_LC(3,:)./ratio_woLC(3,:),[0.5,2.5],3);
grad(4) = unwrap_TDS(freq, ratio_LC(4,:)./ratio_woLC(4,:),[0.5,2.5],3);
delta_n = grad*3e-4/(2*pi*210e-6)

h1=plot(freq, LC_phase(1,:), 'color', colorset(1,:), 'linewidth',2); hold all
h2=plot(freq, LC_phase(2,:), 'color', colorset(2,:), 'linewidth',2);
h3=plot(freq, LC_phase(3,:), 'color', colorset(3,:), 'linewidth',2);
h4=plot(freq, LC_phase(4,:), 'color', colorset(4,:), 'linewidth',2);

xlim([0.2,3.5]); box('on')
xlabel('\omega/2\pi (THz)'); ylabel('\phi_{ps} (бу)')
set(gca, 'fontsize', 20)
legend([h1,h2,h3,h4],{'Pol. 1', 'Pol. 2', ...
    'Pol. 3',  'Pol. 4'},'box', 'off','NumColumns', 2) 

%%
plot_freq = 1.6:0.2:3.4;
plot_ind = interp1(freq, 1:numel(freq), plot_freq, 'nearest');

plot_ratio = ratio_LC(:,plot_ind);
t=linspace(0, 2*pi, 100);
Es = cos(t);

colorset2 = ColorGradient(numel(plot_freq), 'full');
figure('units', 'normalized', 'position', [0.2,0.05,0.5,0.85])
for k=1:4
    subplot(2,2,k)
    for j=1:numel(plot_freq)
        Ep(k,j).E = abs(plot_ratio(k,j)).*cos(t+angle(plot_ratio(k,j)));
        E_abs = sqrt(1+abs(plot_ratio(k,j))^2);
        plot(Ep(k,j).E/E_abs,Es/E_abs, 'linewidth', 2, 'color', [colorset2(j,:),1]); hold all
    end
    S0(k,:) = 1+abs(plot_ratio(k,:)).^2;
    S1(k,:) = abs(plot_ratio(k,:)).^2-1;
    S2(k,:) = 2*abs(plot_ratio(k,:)).*cos(angle(plot_ratio(k,:)));
    S3(k,:) = -2*abs(plot_ratio(k,:)).*sin(angle(plot_ratio(k,:)));
    DOLP(k,:) = sqrt(S1(k,:).^2+S2(k,:).^2)./S0(k,:);
    DOCP(k,:) = S3(k,:)./S0(k,:);
    PCR_L(k,:) = abs(plot_ratio(k,:)/sqrt(2)-1/sqrt(2)).^2./S0(k,:);
    PCR_C1(k,:) = abs(plot_ratio(k,:)/sqrt(2)+1i/sqrt(2)).^2./S0(k,:);
    PCR_C2(k,:) = abs(plot_ratio(k,:)/sqrt(2)-1i/sqrt(2)).^2./S0(k,:);

    xlim([-1,1]); ylim([-1,1]); 
    xlabel('E_p'); ylabel('E_s')
    set(gca,'fontsize', 24);
    daspect([1 1 1])
end

subplot(2,2,1)
ch = colorbar;  colormap(colorset2);
ch.Ticks = [1.6,2.5,3.4];
ch.Location = 'north';
ch.FontSize = 18;
caxis([1.6 3.4]);

%%
figure
subplot(1,2,1)
plot(plot_freq, DOLP(1,:),'o-')
hold all
plot(plot_freq, DOLP(3,:),'^-')
[mean(DOLP(1,:)), mean(DOLP(3,:))]
subplot(1,2,2)
plot(plot_freq, DOCP(2,:),'o-')
hold all
plot(plot_freq, -DOCP(4,:),'^-')
[mean(DOCP(2,:)), mean(DOCP(4,:))]

[mean(PCR_L(1,:)), mean(PCR_C2(2,:)), mean(PCR_C1(4,:))]