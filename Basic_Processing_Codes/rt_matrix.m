function output = rt_matrix(n1,n2,n3,d_array,theta_i,freq,polarization,rt)

% Transfer Function Method based multiple-layer R T model
% Optimized for improving computational speed

if strcmp(rt,'t_multi')
    rt = 't'; % to make it adaptive to nk_solver.m
end

% speed of light
c = 299792458;

if strcmp(rt,'t_single') % must use normal Fresnel's equations to calculate t_single
    theta_i=theta_i*pi/180;
    theta_t12=asin(n1.*sin(theta_i)./n2);
    theta_t23=asin(n1.*sin(theta_i)./n3);
    if polarization=='s'
        t12=2*n1.*cos(theta_i)./(n1.*cos(theta_i)+n2.*cos(theta_t12));
        t23=2*n2.*cos(theta_t12)./(n2.*cos(theta_t12)+n3.*cos(theta_t23));
    elseif polarization=='p'
        t12=2*n1.*cos(theta_i)./(n2.*cos(theta_i)+n1.*cos(theta_t12));
        t23=2*n2.*cos(theta_t12)./(n3.*cos(theta_t12)+n2.*cos(theta_t23));
    end
    beta=2*pi.*freq*10^12.*d_array.*n2.*cos(theta_t12)./c;
    t_single=t12.*t23.*exp(-1i*beta);
    
else % use the transfer matrix method
    % rotate d_array to a column if it's a row
    [m,n] = size(d_array);
    if n > m
        d_array = d_array';
    end

    % number of layers aprat from the first and the last
    n2_layers = numel(n2(:,1));
    d_layers = numel(d_array);
    if n2_layers ~= d_layers
        msg='mismatched number of layers'
        return
    end
    % k vector in free space
    k0 = 2*pi.*freq*10^12./c;
    k0_mat = k0(ones(n2_layers,1),:); % copy to make a matrix

    % coefficient for calculating cosine in each layer
    snell_coef = n1.*sin(theta_i*pi/180);
    % check length and build matrix for r-t coefficients and phase coefficient
    n_length = max([numel(n1), numel(n2(1,:)), numel(n3)]);
    freq_length = numel(freq);
    if n_length > 1 % input n has frequency-dependent values
        % check if it matches the frequency axis
        if n_length ~= freq_length
            msg='mismatched frequency points'
            return
        end
        % extend all refractive indices to be frequency-dependent
        if numel(n1) < n_length
            n1_ext = linspace(n1, n1, n_length);
        else
            n1_ext = n1;
        end
        if numel(n2(1,:)) < n_length
            n2_ext = n2(:, ones(1, n_length));
        else
            n2_ext = n2;
        end
        if numel(n3) < n_length
            n3_ext = linspace(n3, n3, n_length);
        else
            n3_ext = n3;
        end
        % build n matrix, cosine matrix and thickness matrix
        n_matrix = [n1_ext; n2_ext; n3_ext];
        cos_matrix = cos(asin(snell_coef./n_matrix));
        % these 3 must be frequency-dependent for phase coefficient calculation
        n_matrix_full = n_matrix;
        cos_matrix_full = cos_matrix;
        d_matrix = d_array(:, ones(1, n_length));
    else % all n input are frequency-independent
        % n and cosine matrices become an array for r-t coefficient
        n_matrix = [n1; n2; n3];
        cos_matrix = cos(asin(snell_coef./n_matrix));
        % these 3 must be frequency-dependent for phase coefficient calculation
        n_matrix_full = n_matrix(:, ones(1, freq_length));
        cos_matrix_full = cos_matrix(:, ones(1, freq_length));
        d_matrix = d_array(:, ones(1, freq_length));
    end

    % calculate the r-t and phase coefficients as a whole matrix
    n1_used = n_matrix(1:end-1,:);
    n2_used = n_matrix(2:end,:);
    cos1 = cos_matrix(1:end-1,:);
    cos2 = cos_matrix(2:end,:);
    if polarization=='s'
        NC1 = n1_used.*cos1;
        NC2 = n2_used.*cos2;
        r12_matrix = (NC1-NC2)./(NC1+NC2);
        if rt=='t'
            t12_matrix = 2*NC1./(NC1+NC2);
        end
    elseif polarization=='p'
        NC1 = n2_used.*cos1;
        NC2 = n1_used.*cos2;
        r12_matrix = (NC1-NC2)./(NC1+NC2);
        if rt=='t'
            t12_matrix = 2*n1_used.*cos1./(NC1+NC2);
        end
    end

    phase_matrix = k0_mat.*d_matrix.*n_matrix_full(2:end-1,:).*cos_matrix_full(2:end-1,:);
    phase_exp1 = exp(1i*phase_matrix);
    phase_exp2 = exp(-1i*phase_matrix);

    % loop to multiply the transfer matrix of each layer
    MM11 = ones(n2_layers+1, numel(freq)); % initialize MM size
    MM12 = MM11;
    MM21 = MM11;
    MM22 = MM11;
    % transfer matrix for each layer
    if rt == 't'
        M11 = 1./t12_matrix;
        M12 = r12_matrix./t12_matrix;
    else % improve computational speed
        M11 = ones(size(r12_matrix));
        M12 = r12_matrix;
    end
    M21 = M12;
    M22 = M11;

    for k = 1:n2_layers
        if k==1
            MM11(k,:) = M11(k,:);
            MM12(k,:) = M12(k,:);
            MM21(k,:) = M21(k,:);
            MM22(k,:) = M22(k,:);
        else
            MM11(k,:) = MM11(k-1,:).*M11(k,:)+MM12(k-1,:).*M21(k,:);
            MM12(k,:) = MM11(k-1,:).*M12(k,:)+MM12(k-1,:).*M22(k,:);
            MM21(k,:) = MM21(k-1,:).*M11(k,:)+MM22(k-1,:).*M21(k,:);
            MM22(k,:) = MM21(k-1,:).*M12(k,:)+MM22(k-1,:).*M22(k,:);
        end

        P11 = phase_exp1(k,:);
        P22 = phase_exp2(k,:);

        MM11(k,:) = MM11(k,:).*P11;
        MM12(k,:) = MM12(k,:).*P22;
        MM21(k,:) = MM21(k,:).*P11;
        MM22(k,:) = MM22(k,:).*P22;
    end

    % multiply the last transmission matrix, but ONLY MM11 and MM21 will be used
    MM11(k+1,:) = MM11(k,:).*M11(k+1,:)+MM12(k,:).*M21(k+1,:);
end

if rt == 'r'
    MM21(k+1,:) = MM21(k,:).*M11(k+1,:)+MM22(k,:).*M21(k+1,:);
    r_all = MM21(k+1,:)./MM11(k+1,:);
    output = r_all;
elseif rt == 't'
    t_all = 1./MM11(k+1,:);
    output = t_all;
elseif rt == 't_single'
    output = t_single;
end