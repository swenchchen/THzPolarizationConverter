function output=rt_trilayer_solve(n1,n2,n3,n_search,thickness,theta_i,freq,ps,rt)

%%%%%%%%
% calculate the reflection or transmission of a trilayer structure
% for use with 'nk_solution.m'
% contact Swench: swench@qq.com for bugs or problems
%%%%%%%%

% replace the searching nk with the nk matrix
N_size=[numel(n1),numel(n2),numel(n3)];
solve_N_ind=find(N_size==2);
if ~isempty(solve_N_ind)
    if solve_N_ind(1)==1
        n1=n_search;
    elseif solve_N_ind(1)==2
        n2=n_search;
    elseif solve_N_ind(1)==3
        n3=n_search;
    end
end

c=299792458;
theta_i=theta_i*pi/180;
snell_coef = n1.*sin(theta_i);
cos1 = cos(theta_i);
cos2 = cos(asin(snell_coef./n2));
cos3 = cos(asin(snell_coef./n3));

if thickness~=0 % trilayer structure
    if strcmp(ps,'s')
        NC11 = n1.*cos1;
        NC22 = n2.*cos2;
        NC33 = n3.*cos3;
        if strcmp(rt,'r') || strcmp(rt,'t_multi') 
            r12=(NC11-NC22)./(NC11+NC22);
            r23=(NC22-NC33)./(NC22+NC33);
        end
        if strcmp(rt,'t_multi') || strcmp(rt,'t_single')
            t12=2*NC11./(NC11+NC22);
            t23=2*NC22./(NC22+NC33);
        end
    elseif strcmp(ps,'p')
        NC11 = n1.*cos1;
        NC22 = n2.*cos2;
        NC12 = n1.*cos2;
        NC21 = n2.*cos1;
        NC23 = n2.*cos3;
        NC32 = n3.*cos2;
        if strcmp(rt,'r') || strcmp(rt,'t_multi') 
            r12=(NC21-NC12)./(NC21+NC12);
            r23=(NC32-NC23)./(NC32+NC23);
        end
        if strcmp(rt,'t_multi') || strcmp(rt,'t_single')
            t12=2*NC11./(NC21+NC12);
            t23=2*NC22./(NC32+NC23);
        end
    else
        msg='polarization input error'
    end
    
    beta=2*pi.*freq*10^12.*thickness.*NC22./c;
    if strcmp(rt,'r')
        r123=(r12+r23.*exp(-2i*beta))./(1+r12.*r23.*exp(-2i*beta));
        output=r123;
    elseif strcmp(rt,'t_multi')
        t123=(t12.*t23.*exp(-1i*beta))./(1+r12.*r23.*exp(-2i*beta));
        output=t123;
    elseif strcmp(rt,'t_single')
        t123=t12.*t23.*exp(-1i*beta);
        output=t123;
    end
else % bilayer structure, must be reflection (rt='r')
    if strcmp(ps,'s')
        NC11 = n1.*cos1;
        NC33 = n3.*cos3;
        r13=(NC11-NC33)./(NC11+NC33);
    elseif strcmp(ps,'p')
        NC13 = n1.*cos3;
        NC31 = n3.*cos1;
        r13=(NC31-NC13)./(NC31+NC13);
    else
        msg='polarization input error'
    end
    output=r13;
end