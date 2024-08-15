clear all
freq = 0.1:0.01:4;
n1 = 3.419-1i*1e-4;
n2 = 1-1i*1e-2;
metal_sigma = 4*10^7;
d_metal = 200e-9;
metal_prop = nk_epsilon_cond(metal_sigma*d_metal*1000,freq,d_metal,1,'conductivity');
n3 = metal_prop.nk;
thickness = [1,2,3,5,7,10,15,30]*1e-6;
theta_i=30;

figure
subplot(1,2,1)
colorset = ColorGradient(numel(thickness), 'full');
fit_ind = find(freq>1.5 & freq<3.5);
freq_fit = freq(fit_ind);
for k=1:numel(thickness)
    rtp = rt_trilayer(n1,n2,n3,thickness(k),theta_i,freq,'p');
    rp(k,:) = rtp.r;
    rts = rt_trilayer(n1,n2,n3,thickness(k),theta_i,freq,'s');
    rs(k,:) = rts.r;
    phase_diff(k,:) = 180/pi*angle(rp(k,:) ./ rs(k,:));
    phase_fit = phase_diff(k,fit_ind);
    [fitresult, gof] = createFit(freq_fit, phase_fit);
    coeffvals(k,:)= coeffvalues(fitresult);
    rmse_fit(k) = gof.rmse;
    
    plot(freq, 180/pi*angle(rp(k,:) ./ rs(k,:)),'linewidth', 2,'color', colorset(k,:))
    hold all
    plot(freq_fit, coeffvals(k,1)*freq_fit+coeffvals(k,2), 'k--', 'linewidth', 2)
    plot(freq, 180/pi*angle(rp(k,:) ./ rs(k,:))-coeffvals(k,1)*freq, '--','linewidth', 2,'color', [colorset(k,:),0.5])
end
xlim([0.5,4])
set(gca, 'fontsize', 20)
xlabel('f (THz)')
ylabel('Phase diff. (deg.)')
grid('on')
subplot(1,2,2)
yyaxis left
plot(thickness*1e6, coeffvals(:,1), 'o-', 'linewidth',2); hold all
xlabel('Gap thickness (\mum)')
ylabel('p_1')
yyaxis right
plot(thickness*1e6, coeffvals(:,2), '^-', 'linewidth',2)
ylabel('p_2'); set(gca, 'fontsize', 20)

%% 
d_gap = [1.7, 4.6, 12.3, 24]*1e-6;
grad = [38,38,30.7,5];
delta_t = grad/360;

colorset = ColorGradient(7,'7lines');
% figure('units', 'normalized','position',[0.1,0.1,0.65,0.30])
figure('units', 'normalized','position',[0.1,0.1,0.75,0.30])
for k = 1:4
    rtp = rt_trilayer(n1,n2,n3,d_gap(k),theta_i,freq,'p');
    rp(k,:) = rtp.r;
    rts = rt_trilayer(n1,n2,n3,d_gap(k),theta_i,freq,'s');
    rs(k,:) = rts.r;
    phase_diff(k,:) = 180/pi*angle(rp(k,:) ./ rs(k,:));
    phase_LC(k,:) = -freq*grad(k);
    
    subplot(1,4,k)
    h1=plot(freq, phase_diff(k,:), 'k--', 'linewidth', 2); hold all
    h2=plot(freq, phase_LC(k,:), 'k:', 'linewidth', 2);
    phase_diff1 = phase_diff(k,(phase_diff(k,:)<0)==1);  f_phase1 = freq((phase_diff(k,:)<0)==1);
    phase_diff2 = phase_diff(k,(phase_diff(k,:)<0)==0);  f_phase2 = freq((phase_diff(k,:)<0)==0);
    area(f_phase1, phase_diff1,'FaceColor',colorset(1,:),'FaceAlpha',0.1,'linestyle', 'none'); hold all
    area(f_phase2, phase_diff2,'FaceColor',colorset(2,:),'FaceAlpha',0.1,'linestyle', 'none');
    
    phase_LC1 = phase_LC(k,(phase_LC(k,:)<0)==1);  f_phase1 = freq((phase_LC(k,:)<0)==1);
    phase_LC2 = phase_LC(k,(phase_LC(k,:)<0)==0);  f_phase2 = freq((phase_LC(k,:)<0)==0);
    area(f_phase1, phase_LC1,'FaceColor',colorset(1,:),'FaceAlpha',0.11,'linestyle', 'none'); hold all
    area(f_phase2, phase_LC2,'FaceColor',colorset(2,:),'FaceAlpha',0.1,'linestyle', 'none');
    h3=plot(freq, phase_diff(k,:)+phase_LC(k,:), 'linewidth', 4,'color',[0.8,0.1,0.1]);
    xlim([0.8, 3.5]); ylim([-200,120]); yticks([-180,-90,0,90])
    ax = gca;
    ax.YGrid = 'on';
    if k == 1
        xlabel('\omega/2\pi (THz)'); ylabel('\phi_{ps} (бу)');
    else
        yticklabels([]);
    end
    set(gca, 'fontsize', 18); 
end
legend([h1, h2, h3], {'mirror-TIR', 'LC', 'Combined'}, 'fontsize', 16,'location','south','box','off');


%%
function [fitresult, gof] = createFit(freq_fit, phase_fit)

[xData, yData] = prepareCurveData( freq_fit, phase_fit );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );
end