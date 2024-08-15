function output=nk_solution(freq,ratio,Ns1,Ns2,Ns3,Nr1,Nr2,Nr3,ds2,dr2,rt_s,rt_r,theta_i,ps_s,ps_r,accuracy)
%%  INSTRUCTIONS
% ----------------------------------------------LAST UPDATE: 2022.05.09 ------------------------------------------------------------
%  This function solves 1 set of unknown nk for most general THz measurements
%  Including transmission (single transmission & mutiple transmission), reflection (normal reflection and ellipsometry).
%  Please contact Swench (swench@qq.com) for any problem with the code
%
% [freq]: only the frequency region of interests, in THz unit, eg. 0.1-2 THz
%
% [ratio]: ratio=Esample/Ereference in frequency-domain (experimental data)
%
% [Ns1, Ns2, Ns3]: the complex refractive indices of the trilayer optical model for measuring Esample. 
% Input as a value if frequency-independent
% Input as a 1D array if frequency-dependent. 
% Input as a range for the FIRST FREQUENCY POINT if unknown. Set a LARGER range if unsure e.g. Ns2=[1-1i*0,10-1i*10]
%
% [Nr1, Nr2, Nr3]: the complex refractive indices of the trilayer optical model for measuring Ereference.
% May also contain unknown variable in ellipsometry measurement
%
% [ds2,dr2]: thickness of the 2nd layer for sample (ds2) and reference (dr2) respectively, in unit of m
%
% [rt_s,rt_r]: measurement mode for sample (rt_s) and reference (rt_r)
% ='r': reflection
% ='t_single': only consider the 1st transmission.  
%    If both rt_s and rt_r are 't_single', Beer-Lambert's Law (analytical solution) will be used. 
%    In this case, the program regards dr2=ds2, Nr1=Ns1, Nr2=1, Nr3=Ns3,  theta_i=0, accuracy: not applicable
% ='t_multi': consider multiple transmission
%
% [theta_i]: incident angle in degree
%
% [ps_s] [ps_r] = 'p' or 's', p- or s- polarizations for sample (ps1) and reference (ps2)
%
% [accuracy] how small the relative error between the calculated ratio and the measured ratio. Suggesting value: 0.01~0.001

%% Recognize the unknwon nk to be solved
N_size=[numel(Ns1),numel(Ns2),numel(Ns3),numel(Nr1),numel(Nr2),numel(Nr3)];
solve_N_ind=find(N_size==2);

if solve_N_ind(1)==1;   solve_N=Ns1;
elseif solve_N_ind(1)==2;   solve_N=Ns2;
elseif solve_N_ind(1)==3;   solve_N=Ns3;
elseif solve_N_ind(1)==4;   solve_N=Nr1;
elseif solve_N_ind(1)==5;   solve_N=Nr2;
elseif solve_N_ind(1)==6;   solve_N=Nr3;
end

if max(solve_N_ind)>3 % Nr contains unknown values to be solved
    unchanged_ref = 0; % E_ref (E2) should be calculated in every searching loop
else
    unchanged_ref = 1;
end
%% Solution start
c = 299792458; % speed of light in m/s
alpha_k_ratio=400*pi.*freq./(c*10^-8);

if strcmp(rt_s,'t_single') && strcmp(rt_r,'t_single') % analytical solution by Beer-Lambert's Law
    phase_unwrap_raw=unwrap(angle(ratio));
    poly_coef = polyfit(freq, phase_unwrap_raw, 1);
    poly_shift = poly_coef(2); %take absolute shift
    num_pi = round(poly_shift ./ (2.*pi)); 
    phase_shift = num_pi .* 2 .* pi;
    phase_unwrap=phase_unwrap_raw-phase_shift;
    
    n_output=1 - phase_unwrap.*c./(2*pi*freq*10^12.*ds2);
    t_coef1=2*Ns1./(Ns1+n_output);
    t_coef2=2*n_output./(Ns3+n_output);
    t_coef=t_coef1.*t_coef2;
    a_output=-0.02*log(abs(ratio./t_coef)) / ds2; % alpha in cm-1
    k_output=a_output./alpha_k_ratio;
    
else % numerical solution for other optical models
    mesh_num_init=100; % a large nk matrix 200*200 for the first frequency point to ensure the convergency
    mesh_num=20;
    n_range_initial=linspace(real(solve_N(1)),real(solve_N(2)),mesh_num_init);
    k_range_initial=linspace(-imag(solve_N(1)),-imag(solve_N(2)),mesh_num_init);

    n_output=linspace(0,0,numel(freq));
    k_output=linspace(0,0,numel(freq));

    freq_step=freq(2)-freq(1);
    failed_flag=0;
    
    if unchanged_ref == 1
        E2_full=rt_trilayer_solve(Nr1,Nr2,Nr3,0,dr2,theta_i,freq,ps_r,rt_r);
    end
    
    for k=1:numel(freq)
        current_f=freq(k);
        current_ratio=ratio(k);
        % if N size not equal to 2, it's input as a known parameter,
        % min(...) extracts the current frequency point, or the first point if it's frequency-independent.
        if numel(Ns1)~=2;  Ns1_used=Ns1(min(numel(Ns1),k));  else;  Ns1_used=Ns1;  end
        if numel(Ns2)~=2;  Ns2_used=Ns2(min(numel(Ns2),k));  else;  Ns2_used=Ns2;  end
        if numel(Ns3)~=2;  Ns3_used=Ns3(min(numel(Ns3),k));  else;  Ns3_used=Ns3;  end
        if unchanged_ref == 1
            E2 = E2_full(min([numel(E2_full),k]));
        else
            if numel(Nr1)~=2;  Nr1_used=Nr1(min(numel(Nr1),k));  else;  Nr1_used=Nr1;  end
            if numel(Nr2)~=2;  Nr2_used=Nr2(min(numel(Nr2),k));  else;  Nr2_used=Nr2;  end
            if numel(Nr3)~=2;  Nr3_used=Nr3(min(numel(Nr3),k));  else;  Nr3_used=Nr3;  end
        end
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
            if n_ref<10  % it should be a dielectric medium
                nk_step=8*freq_step; % the nk searching range depends on the frequency step. 8*freq_step is according to the dispersion of water
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
            E1=rt_trilayer_solve(Ns1_used,Ns2_used,Ns3_used,nk_matrix,ds2,theta_i,current_f,ps_s,rt_s);
            if unchanged_ref == 0
                E2=rt_trilayer_solve(Nr1_used,Nr2_used,Nr3_used,nk_matrix,dr2,theta_i,current_f,ps_r,rt_r);
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