function output = rt_matrix_single(n_all, solve_num, nk_matrix, d_array,theta_i,freq,polarization,rt)

%%%%%%%%%%%% Principles %%%%%%%%%%%%
% RT_out = M1*M2*M3*M4*M5*M6*M7

% a loop to calculate the above multiplication

% the difference compared to rt_matrix, is each M is individually
% calculated in each loop here, while in rt_matrix calculats all the
% relavant matrices first, and only read the values in each loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(rt,'t_multi')
    rt = 't'; % to make it adaptive to nk_solver.m
end

% speed of light
c = 299792458;

% rotate d_array to a column if it's a row
[m,n] = size(d_array);
if n > m
    d_array = d_array';
end

% number of layers aprat from the first and the last
n_mid_layers = numel(n_all)-2;
d_layers = numel(d_array);
if n_mid_layers ~= d_layers
    msg='mismatched number of layers'
    return
end

% k vector in free space, single value
k0 = 2*pi.*freq*10^12./c;

% cosine in each layer
snell_coef = n_all(1).*sin(theta_i*pi/180);
cos_all = cos(asin(snell_coef./n_all));
cos_matrix = cos(asin(snell_coef./nk_matrix));

for k = 1:numel(d_array)+1
    if k~= solve_num-1 && k~=solve_num % sample not involved
        n_left = n_all(k);
        n_right = n_all(k+1);
        cos_left = cos_all(k);
        cos_right = cos_all(k+1);
    elseif k== solve_num-1 % sample on the right
        n_left = n_all(k);
        n_right = nk_matrix;
        cos_left = cos_all(k);
        cos_right = cos_matrix;
    elseif k== solve_num % sample on the left
        n_left = nk_matrix;
        n_right = n_all(k+1);
        cos_left = cos_matrix;
        cos_right = cos_all(k+1);
    end
        
    if polarization=='s'
        NC1 = n_left.*cos_left;
        NC2 = n_right.*cos_right;
        r12_matrix = (NC1-NC2)./(NC1+NC2);
        if rt=='t'
            t12_matrix = 2*NC1./(NC1+NC2);
        end
    elseif polarization=='p'
        NC1 = n_right.*cos_left;
        NC2 = n_left.*cos_right;
        r12_matrix = (NC1-NC2)./(NC1+NC2);
        if rt=='t'
            t12_matrix = 2*n_left.*cos_left./(NC1+NC2);
        end
    end
    
    % phase
    if k~=numel(d_array)+1 % i.e. no phase matrix for the last layer 
        phase_matrix = k0.*d_array(k).*n_right.*cos_right;
        P11 = exp(1i*phase_matrix);
        P22 = exp(-1i*phase_matrix);
    end
    
    if rt == 't'
        M11 = 1./t12_matrix;
        M12 = r12_matrix./t12_matrix;
    else % just assume t12=1 to improve computational speed
        M11 = ones(size(r12_matrix));
        M12 = r12_matrix;
    end
    M21 = M12;
    M22 = M11;
    
    if k==1
        MM11 = M11;
        MM12 = M12;
        MM21 = M21;
        MM22 = M22;
    elseif k~=numel(d_array)+1
        MM11_new = MM11.*M11+MM12.*M21;
        MM12_new = MM11.*M12+MM12.*M22;
        MM21_new = MM21.*M11+MM22.*M21;
        MM22_new = MM21.*M12+MM22.*M22;
        
        MM11 = MM11_new;
        MM12 = MM12_new;
        MM21 = MM21_new;
        MM22 = MM22_new;
    else % the last layer, only need MM11 & MM21 (reflection only)
        MM11 = MM11.*M11+MM12.*M21;
        if rt == 'r'
            MM21 = MM21.*M11+MM22.*M21;
        end
    end
    
    % transfer matrix * phase matrix
    if k~=numel(d_array)+1
        MM11 = MM11.*P11;
        MM12 = MM12.*P22;
        MM21 = MM21.*P11;
        MM22 = MM22.*P22;
    end
end

if rt == 'r'
    r_all = MM21./MM11;
    output = r_all;
elseif rt == 't'
    t_all = 1./MM11;
    output = t_all;
end