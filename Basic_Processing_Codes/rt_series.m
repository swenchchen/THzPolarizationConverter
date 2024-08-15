function output=rt_series(n1,n_matrix,d_array,theta_i,freq,polarization,rt)
% n_matrix should be >=3 lines
% d_array should be >=2 values,corresponding to d of n_matrix(1:end-1,:)
% the sequence is from n1 to deeper layers

%%
c=299792458;
theta_i=theta_i*pi/180;
n_sin_product=n1.*sin(theta_i);
theta_t_matrix=asin(n_sin_product./n_matrix);

r_fmb=zeros(numel(d_array),numel(freq));
t_fmb=zeros(numel(d_array),numel(freq));
for k=1:numel(d_array)-1
    
    n_back=n_matrix(end-k+1,:);
    n_middle=n_matrix(end-k,:);
    n_former=n_matrix(end-k-1,:);
    
    theta_t_back=theta_t_matrix(end-k+1,:);
    theta_t_middle=theta_t_matrix(end-k,:);
    theta_t_former=theta_t_matrix(end-k-1,:);
    
    d=d_array(end-k+1);
    
    if k==1
        if polarization=='s'
            r_fm=(n_former.*cos(theta_t_former)-n_middle.*cos(theta_t_middle))./(n_former.*cos(theta_t_former)+n_middle.*cos(theta_t_middle));
            r_mb=(n_middle.*cos(theta_t_middle)-n_back.*cos(theta_t_back))./(n_middle.*cos(theta_t_middle)+n_back.*cos(theta_t_back));
            if rt=='t'
                t_fm=2*n_former.*cos(theta_t_former)./(n_former.*cos(theta_t_former)+n_middle.*cos(theta_t_middle));
                t_mb=2*n_middle.*cos(theta_t_middle)./(n_middle.*cos(theta_t_middle)+n_back.*cos(theta_t_back));
            end
        elseif polarization=='p'
            r_fm=(n_middle.*cos(theta_t_former)-n_former.*cos(theta_t_middle))./(n_middle.*cos(theta_t_former)+n_former.*cos(theta_t_middle));
            r_mb=(n_back.*cos(theta_t_middle)-n_middle.*cos(theta_t_back))./(n_back.*cos(theta_t_middle)+n_middle.*cos(theta_t_back));
            if rt=='t'
                t_fm=2*n_former.*cos(theta_t_former)./(n_middle.*cos(theta_t_former)+n_former.*cos(theta_t_middle));
                t_mb=2*n_middle.*cos(theta_t_middle)./(n_back.*cos(theta_t_middle)+n_middle.*cos(theta_t_back));
            end
        end
        beta_middle=2*pi.*freq*10^12.*d.*n_middle.*cos(theta_t_middle)./c;

        r_fmb(k,:)=(r_fm+r_mb.*exp(-2i*beta_middle))./(1+r_fm.*r_mb.*exp(-2i*beta_middle));
        if rt=='t'
            t_fmb(k,:)=(t_fm.*t_mb.*exp(-1i*beta_middle))./(1+r_fm.*r_mb.*exp(-2i*beta_middle));
        end
    else
        if polarization=='s'
            r_fm=(n_former.*cos(theta_t_former)-n_middle.*cos(theta_t_middle))./(n_former.*cos(theta_t_former)+n_middle.*cos(theta_t_middle));
            if rt=='t'
                t_fm=2*n_former.*cos(theta_t_former)./(n_former.*cos(theta_t_former)+n_middle.*cos(theta_t_middle));
            end
        elseif polarization=='p'
            r_fm=(n_middle.*cos(theta_t_former)-n_former.*cos(theta_t_middle))./(n_middle.*cos(theta_t_former)+n_former.*cos(theta_t_middle));
            if rt=='t'
                t_fm=2*n_former.*cos(theta_t_former)./(n_middle.*cos(theta_t_former)+n_former.*cos(theta_t_middle));
            end
        end
        beta_middle=2*pi.*freq*10^12.*d.*n_middle.*cos(theta_t_middle)./c;
        
        r_fmb(k,:)=(r_fm+r_fmb(k-1,:).*exp(-2i*beta_middle))./(1+r_fm.*r_fmb(k-1,:).*exp(-2i*beta_middle));
        if rt=='t'
            t_fmb(k,:)=(t_fm.*t_fmb(k-1,:).*exp(-1i*beta_middle))./(1+r_fm.*r_fmb(k-1,:).*exp(-2i*beta_middle));
        end
    end
end

% last trilayer including the window/substract
n_middle=n_matrix(end-k-1,:);
n_former=n1;
theta_t_middle=theta_t_matrix(end-k-1,:);
theta_t_former=theta_i;
d=d_array(end-k);

if polarization=='s'
    r_fm=(n_former.*cos(theta_t_former)-n_middle.*cos(theta_t_middle))./(n_former.*cos(theta_t_former)+n_middle.*cos(theta_t_middle));
    if rt=='t'
        t_fm=2*n_former.*cos(theta_t_former)./(n_former.*cos(theta_t_former)+n_middle.*cos(theta_t_middle));
    end
elseif polarization=='p'
    r_fm=(n_middle.*cos(theta_t_former)-n_former.*cos(theta_t_middle))./(n_middle.*cos(theta_t_former)+n_former.*cos(theta_t_middle));
    if rt=='t'
        t_fm=2*n_former.*cos(theta_t_former)./(n_middle.*cos(theta_t_former)+n_former.*cos(theta_t_middle));
    end
end
beta_middle=2*pi.*freq*10^12.*d.*n_middle.*cos(theta_t_middle)./c;

r_fmb(k+1,:)=(r_fm+r_fmb(k,:).*exp(-2i*beta_middle))./(1+r_fm.*r_fmb(k,:).*exp(-2i*beta_middle));
if rt=='t'
    t_fmb(k+1,:)=(t_fm.*t_fmb(k,:).*exp(-1i*beta_middle))./(1+r_fm.*r_fmb(k,:).*exp(-2i*beta_middle));
end

if rt=='r'
    output=r_fmb(k+1,:);
else
    output=t_fmb(k+1,:);
end

