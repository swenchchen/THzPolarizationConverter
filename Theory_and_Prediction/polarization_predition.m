clear all

% Input the target polarization
theta = 45;
epsilon = -0.01;
target_freq = 2;
% save_file = 'phase_15_n40.mat';
tolerance = 8;

% Convert to Stokes Vector
S1 = cos(2*epsilon*pi/180)*cos(2*theta*pi/180);
S2 = cos(2*epsilon*pi/180)*sin(2*theta*pi/180);
S3 = sin(2*epsilon*pi/180);
sin2Psi = sqrt(1-S1^2);

Delta1 = [acos(S2/sqrt(1-S1^2)), -acos(S2/sqrt(1-S1^2))]; 
S3_potential = -sin2Psi*sin(Delta1); 
min_ind = find(abs(S3_potential-S3) == min(abs(S3_potential-S3))); 
Delta = Delta1(min_ind(1)); 
Psi = 0.5*acos(-S1); 
tanPsi = tan(Psi);
target_phase = Delta*180/pi; 

%%
if target_phase >95 
    % Cannot reach >95, need -180
    target_phase = target_phase-180;
    input_pol = 180;
elseif target_phase >0
    % two options
    target_phase = [target_phase, target_phase-180];
    input_pol = [0, 180];
elseif target_phase >-90
    % one option, cannot +180
    input_pol = 0;
else
    % two options
    target_phase = [target_phase, target_phase+180];
    input_pol = [0, 180];
end

freq = 0.1:0.01:20; % searching range
freq_step = freq(2)-freq(1);
target_f_ind = find(abs(freq-target_freq) == min(abs(freq-target_freq)));

n1 = 3.418-1i*1e-4;
n2 = 1-1i*1e-2;
metal_sigma = 4*10^7;
d_metal = 200e-9;
metal_prop = nk_epsilon_cond(metal_sigma*d_metal*1000,freq,d_metal,1,'conductivity');
n3 = metal_prop.nk;
theta_i=30;

loop_num = 3;
for t = 1:numel(target_phase) 
    clear match_num
    for j = 1:loop_num
        % 3 loops to find the required height and Delta_t
        if j==1
            % rough search in the first loop
            thickness = 0.5*1.3.^(1:20)*1e-6;
        elseif j == 2
            % finner scan according to the result from the first loop
            thickness = linspace(0.7*best_thickness, 1.5*best_thickness, 20);
        elseif j==3
            % calculate multiple combinations of height values & Delta t,
            % compare the phase and get the optimized solution
            thickness = linspace(best_thickness, 4*best_thickness, 200);
        end

        phase_comp = zeros(numel(thickness), numel(freq));
        phase_diff = zeros(numel(thickness), numel(freq));
        for k=1:numel(thickness)
            rtp = rt_trilayer(n1,n2,n3,thickness(k),theta_i,freq,'p');
            rts = rt_trilayer(n1,n2,n3,thickness(k),theta_i,freq,'s');
            phase_diff(k,:) = 180/pi*angle(rtp.r ./ rts.r);
            phase_grad = diff(phase_diff(k,:));
            % Gradient and the center of the frequency range
            phase_grad_center = phase_grad(target_f_ind)/freq_step;
            if j<3
                % Roughly calculate the required phase gradient：freq*phase_grad_center
                phase_comp(k,:) = phase_diff(k,:) - freq*phase_grad_center;
            else
                % The third loop，need to finely scan different Delta_t
                grad_var = 1.5*tolerance/target_freq;
                % Gradient range
                grad_range = linspace(phase_grad_center-0.5*grad_var, phase_grad_center+2*grad_var, 100);
                % Loop for every gradient range
                for m = 1:numel(grad_range)
                    phase_comp_grad = phase_diff(k,:) - freq*grad_range(m);
                    phase_err = abs(phase_comp_grad-target_phase(t));
                    % Within the tolerance
                    match_freq = freq(phase_err<=tolerance);
                    if numel(match_freq)>1 % there's at least one solution
                        if max(diff(match_freq))>2*freq_step || min(abs(match_freq-target_freq))>freq_step
                            % inconsistent range, or the target frequency
                            % falls out of the range
                            match_freq = []; % returns zero matching frequency
                        end
                    end

                    if isempty(match_freq)
                        match_range = 0;
                    else
                        if match_freq(1)>freq(1) && match_freq(end)<freq(end)
                            % The range may not centered at the required
                            % center
                            % So calculate the side having the shortest
                            % distance to the required center
                            match_range = min([target_freq - match_freq(1), match_freq(end)-target_freq]);
                        else
                            % the range is very hugh, will take the whole
                            % range into consideration, ignore the center
                            % position
                            match_range = match_freq(end)-match_freq(1);
                        end
                    end

                    match_num(k,m) = match_range;
                end
            end
        end
        thickness_ind = find(abs(phase_comp(:,target_f_ind)-target_phase(t)) == min(abs(phase_comp(:,target_f_ind)-target_phase(t))));
        best_thickness = thickness(thickness_ind);
    end

%     figure
%     imagesc(match_num)
    [k_best, m_best] = find(match_num == max(match_num,[],'all'));

    phase_diff_best = phase_diff(k_best(1),:);
    phase_grad = diff(phase_diff_best);
    phase_grad_center = phase_grad(target_f_ind)/freq_step;
    grad_range = linspace(phase_grad_center-0.5*grad_var, phase_grad_center+2*grad_var, 100); %要保持和上面更改同步！！
    thickness_best(t) = thickness(k_best(1));
    delta_t_best(t) = grad_range(m_best(1))/360;

    phase_comp_best(t,:) = phase_diff_best - freq*grad_range(m_best(1));
    phase_err = abs(phase_comp_best(t,:)-target_phase(t));
    match_ind = find(phase_err<tolerance);
    match_bandwidth(t,:) = [freq(match_ind(1)), freq(match_ind(end))];
    
    if t==1
        figure
    end
%     plot(freq, phase_comp_best(t,:))
%     hold all
    plot(freq(match_ind), phase_comp_best(t,match_ind), 'o'); hold all
    ylim([-200,130])
%     title(['d=',num2str(thickness_best), ', t=',num2str(delta_t_best)]);
    phase_predict(t).data(1,:) = freq(match_ind);
    phase_predict(t).data(2,:) = phase_comp_best(t,match_ind);
end
% save(save_file, 'phase_predict');

['Ep/Es=', num2str(real(exp(1i*input_pol*pi/180).*tanPsi))]
['gap thickness = ', num2str(thickness_best*1e6), 'um']
['t = ',num2str(delta_t_best), 'ps']
if numel(match_bandwidth(:,1))<2 % Only one option
    ['bandwidth = ', num2str(match_bandwidth)]
else
    ['bandwidth = ', num2str([match_bandwidth(1,:),match_bandwidth(2,:)])]
end