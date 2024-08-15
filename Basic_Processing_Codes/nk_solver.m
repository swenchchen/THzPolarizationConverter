function output = nk_solver(freq, ratio, N_sam, N_ref, solve_num, d_sam, d_ref, theta_i, rt2, ps2, accuracy)
%%  INSTRUCTIONS
% 
%  Upgraded version of nk_solution, being able to solve optical models with
%  arbitrary number of layers
%
%  Please contact Swench (swench@qq.com) for any problem about the code
%
% [freq]: only the frequency region of interests, in THz unit, eg. 0.1-2 THz
%
% [ratio]: ratio=Esample/Ereference in frequency-domain (experimental data)
%
% N_sam is a array structure, N_sam(x).n is the complex refractive index of x layer
% Same for N_ref
% 
% solve_num is an array of 2 values [a,b], indicating the unknown nk layer
% number in N_sam and N_ref
%
% d_sam is an array of thicknesses for all layers with finite thickness,
% Make sure numel(d_sampl) = numel(N_sam)-2
% Same for d_ref
%
% [theta_i]: incident angle in degree
%
% rt2 is an array of structure: rt2(1).rt and rt2(2).rt for measurement modes of sample and reference
% ='r': reflection
% ='t_single': only consider the 1st transmission.  
%    If both rt_s and rt_r are 't_single', Beer-Lambert's Law (analytical solution) will be used. 
%    In this case, the program regards dr2=ds2, Nr1=Ns1, Nr2=1, Nr3=Ns3,  theta_i=0, accuracy: not applicable
% ='t_multi': consider multiple transmission
%
%  ps2 is an array of structure: ps2(1).ps and ps2(2).ps for polarizations of sample and reference
% = 'p' or 's'
%
% [accuracy] how small the relative error between the calculated ratio and
% the measured ratio. Recommended value: 1e-2~1e-4;
%
% Template:
% N_sam(1).n = 1-1i*1e-4;  N_sam(2).n = 1-1i*1e-4;  N_sam(3).n = 1-1i*1e-4;  
% N_ref(1).n = 1-1i*1e-4;  N_ref(2).n = 1-1i*1e-4;  N_ref(3).n = 1-1i*1e-4;  
% solve_num = [2,0];  d_sam = 5e-6; d_ref = d_sam;  theta_i = 0;
% rt2(1).rt = 'r';  rt2(2).rt = 'r';  ps2(1).ps = 'p';  ps2(2).ps = 'p';  accuracy = 1e-3;

%% Recognize the unknwon nk to be solved
solve_N = N_sam(solve_num(1)).n;

if solve_num(2) ~=0 % N_ref contains unknown values to be solved
    unchanged_ref = 0; % E_ref (E2) should be calculated in every searching loop
else
    unchanged_ref = 1;
end

%% Solution start
c = 299792458; % speed of light in m/s
alpha_k_ratio=400*pi.*freq./(c*10^-8);

if strcmp(rt2(1).rt,'t_single') && strcmp(rt2(2).rt,'t_single') % analytical solution by Beer-Lambert's Law
    % Both sample and reference must be tri-layer structures, with the N_sam(2) as the unknown
    fit_range = find(freq >0.5 & freq<1.5); % disable this if want to unwrap according to the full range of [freq]
%     fit_range = 1:numel(freq);% enalbe this if want to unwrap according to the full range of [freq]
    phase_unwrap_raw=unwrap(angle(ratio));
    poly_coef = polyfit(freq(fit_range), phase_unwrap_raw(fit_range), 1);
    poly_shift = poly_coef(2); %take absolute shift
    num_pi = round(poly_shift ./ (2.*pi)); 
    phase_shift = num_pi .* 2 .* pi;
    phase_unwrap = phase_unwrap_raw-phase_shift;
    
    n_output = real(N_ref(2).n) - phase_unwrap.*c./(2*pi*freq*10^12.*d_sam);
    
    t_coef_sam1=2*N_sam(1).n./(N_sam(1).n+n_output);
    t_coef_sam2=2*n_output./(N_sam(3).n+n_output);
    t_coef_sam=t_coef_sam1.*t_coef_sam2;
    
    t_coef_ref1=2*N_ref(1).n./(N_ref(1).n+N_ref(2).n);
    t_coef_ref2=2*N_ref(2).n./(N_sam(3).n+N_ref(2).n);
    t_coef_ref=t_coef_ref1.*t_coef_ref2;
    
    a_output=-0.02*log(abs(ratio.*t_coef_ref./t_coef_sam)) / d_sam; % alpha in cm-1
    k_output=a_output./alpha_k_ratio;
    
else % numerical solution for other optical models
    
    % extend all layers to be frequency-dependent
    freq_length = numel(freq);
    % reference N
    N_ref_mid = zeros(numel(N_ref)-2, numel(freq));
    N_ref_all = zeros(numel(N_ref),numel(freq));
    for k=1:numel(N_ref)
        if numel(N_ref(k).n) < freq_length
            N_ref(k).n = ones(1,freq_length)*N_ref(k).n;
        end
        if solve_num(2) == 0
            if k>1 && k<numel(N_ref)
                N_ref_mid(k-1,:) = N_ref(k).n;
            end
        else
            N_ref_all(k,:) = N_ref(k).n;
        end
    end
    % sample N
    N_sam_all = zeros(numel(N_sam), numel(freq));
    for k=1:numel(N_sam)
        if numel(N_sam(k).n) < freq_length
            N_sam(k).n = ones(1,freq_length)*N_sam(k).n;
        end
        N_sam_all(k,:) = N_sam(k).n;
    end

    if solve_num(2) == 0 % E_ref is unchanged, only needs to be calculated once
        E2_full=rt_matrix(N_ref(1).n, N_ref_mid, N_ref(end).n, d_ref, theta_i, freq, ps2(2).ps, rt2(2).rt);
    end
    
    mesh_num_init = 100; % a large nk matrix 200*200 for the first frequency point to ensure the convergency
    mesh_num = 20;
    
    if real(solve_N)<8 % n is small, set +-3 range for Re(n) and Im(n)
        n_range_initial = linspace(max([real(solve_N)-0.5,1]),  real(solve_N)+0.5,mesh_num_init);
        k_range_initial=linspace(max([-imag(solve_N)-0.1, -0.2]),  -imag(solve_N)+0.1,mesh_num_init);
    else % n is large, set a relative range for Re(n) and Im(n)
        n_range_initial = linspace(real(solve_N)*0.5,  real(solve_N)*2,mesh_num_init);
        k_range_initial=linspace(-imag(solve_N)*0.5,  -imag(solve_N)*2,mesh_num_init);
    end

    n_output=linspace(0,0,numel(freq));
    k_output=linspace(0,0,numel(freq));

    freq_step=freq(2)-freq(1);
    failed_flag=0;
    
    for k=1:numel(freq)

        current_f=freq(k);
        current_ratio=ratio(k);
        
        % assign the searching range for the current frequency
        if k==1 || k-1-failed_flag<1  % It's the 1st frequency point, or it has never successfully found the first solution
            [n_matrix,k_matrix]=meshgrid(n_range_initial,k_range_initial);
            % search in the initial range with a fine mesh
            nk_matrix=n_matrix-1i*k_matrix;
            n_range=n_range_initial;
            k_range=k_range_initial;
        else % search in a nk range according to the value of the last successful frequency point
            n_ref=n_output(k-1-failed_flag);
            k_ref=k_output(k-1-failed_flag);
            if n_ref<8  % it should be a dielectric medium
                if strcmp(rt2(1).rt,'t_multi') && d_sam>100e-6 && d_sam<500e-6
                    nk_step=8*freq_step/(d_sam/100e-6); % reduce the step according to the thickness to avoid local minimums
                elseif strcmp(rt2(1).rt,'t_multi') && d_sam>=500e-6
                    nk_step=8*freq_step/5; % if it's too thick, set the step to a constant value
                else
                    nk_step=8*freq_step; % the nk searching range depends on the frequency step. 8*freq_step is according to the dispersion of water
                end
                n_range=linspace(n_ref-nk_step,  n_ref+nk_step,  mesh_num);
                k_range=linspace(k_ref-nk_step,  k_ref+nk_step,  mesh_num);
            else  % it's possibly a conductive material, change the absolute step to a magnification factor of 1.2 and 0.8
                n_range=linspace(max(n_ref*0.8, 0.8),  n_ref*1.2,  mesh_num);
                if abs(k_ref>3)
                    k_range=linspace(max(k_ref*0.8,-200/418/current_f),  k_ref*1.2,  mesh_num);
                else %if k is too small the solution may exist at k<0, need -step to get into the negative domain
                    nk_step=8*freq_step;
                    k_range=linspace(k_ref-nk_step,  k_ref+nk_step,  mesh_num);
                end
            end
            [n_matrix,k_matrix]=meshgrid(n_range,k_range);
            nk_matrix=n_matrix-1i*k_matrix;
        end
        
        % seaching loop start
        max_searching_num=10;
        for search_loop=1:max_searching_num % maximum 10 times zooming in searching space
            E1 =  rt_matrix_single(N_sam_all(:,k), solve_num(1), nk_matrix, d_sam, theta_i, current_f, ps2(1).ps, rt2(1).rt);
            if solve_num(2) ~= 0 % E2 contains unknow values
                E2 =  rt_matrix_single(N_ref_all(:,k), solve_num(2), nk_matrix, d_ref, theta_i, current_f, ps2(2).ps, rt2(2).rt);
            else
                E2 = E2_full(k);
            end
            ratio_cal=E1./E2;
            ratio_diff=abs(ratio_cal-current_ratio);
            min_ratio_diff=min(ratio_diff,[],'all');
            [k_min_ind,n_min_ind]=find(ratio_diff==min_ratio_diff,1);
            ratio_accuracy=min_ratio_diff/abs(current_ratio);
            if ratio_accuracy>accuracy && search_loop<max_searching_num % haven't reached accuracy requirement, zoom in nk space
                if k_min_ind==1 || k_min_ind==numel(k_range) || n_min_ind==1 || n_min_ind==numel(n_range)
                    % solution found at the boundary, move the searching region without zooming in
                    k_length = (k_range(end) - k_range(1));
                    n_length = (n_range(end) - n_range(1));
                    k_range=linspace(k_range(k_min_ind)-k_length,k_range(k_min_ind)+k_length,numel(k_range));
                    % use numel(k_range) rather than mesh_num because it could be the first point using mesh_num_initial
                    n_range=linspace(n_range(n_min_ind)-n_length,n_range(n_min_ind)+n_length,numel(n_range));
                else % zoom in only when both n and k locate within the searching space
                    k_step=k_range(2)-k_range(1);
                    n_step=n_range(2)-n_range(1);
                    k_range=linspace(k_range(k_min_ind)-2*k_step,k_range(k_min_ind)+2*k_step,mesh_num);
                    n_range=linspace(n_range(n_min_ind)-2*n_step,n_range(n_min_ind)+2*n_step,mesh_num);
                end
                [n_matrix,k_matrix]=meshgrid(n_range,k_range);
                nk_matrix=n_matrix-1i*k_matrix;
            elseif ratio_accuracy<=accuracy
                n_output(k)=n_range(n_min_ind);
                k_output(k)=k_range(k_min_ind);
                failed_flag=0;
                break % stop searching
            else
                n_output(k)=n_range(n_min_ind);
                k_output(k)=k_range(k_min_ind);
                failed_flag=failed_flag+1;
                msg='Failed to optimize within 10 iteractions'
            end
        end
    end
    a_output=k_output.*alpha_k_ratio;
end


%% output results
output.n=n_output;
output.k=k_output;
output.a=a_output;