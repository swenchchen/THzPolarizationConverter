function output = rt_trilayer(n1,n2,n3,thickness,theta_i,freq,polarization)
%%
c=299792458;
theta_i=theta_i*pi/180;
theta_t12=asin(n1.*sin(theta_i)./n2);
theta_t23=asin(n1.*sin(theta_i)./n3);

if polarization=='s'
    r12=(n1.*cos(theta_i)-n2.*cos(theta_t12))./(n1.*cos(theta_i)+n2.*cos(theta_t12));
    r23=(n2.*cos(theta_t12)-n3.*cos(theta_t23))./(n2.*cos(theta_t12)+n3.*cos(theta_t23));
    t12=2*n1.*cos(theta_i)./(n1.*cos(theta_i)+n2.*cos(theta_t12));
    t23=2*n2.*cos(theta_t12)./(n2.*cos(theta_t12)+n3.*cos(theta_t23));
elseif polarization=='p'
    r12=(n2.*cos(theta_i)-n1.*cos(theta_t12))./(n2.*cos(theta_i)+n1.*cos(theta_t12));
    r23=(n3.*cos(theta_t12)-n2.*cos(theta_t23))./(n3.*cos(theta_t12)+n2.*cos(theta_t23));
    t12=2*n1.*cos(theta_i)./(n2.*cos(theta_i)+n1.*cos(theta_t12));
    t23=2*n2.*cos(theta_t12)./(n3.*cos(theta_t12)+n2.*cos(theta_t23));
else
    msg='polarization input error'
end

beta=2*pi.*freq*10^12.*thickness.*n2.*cos(theta_t12)./c;

r123=(r12+r23.*exp(-2i*beta))./(1+r12.*r23.*exp(-2i*beta));
t123=(t12.*t23.*exp(-1i*beta))./(1+r12.*r23.*exp(-2i*beta));
t_single=t12.*t23.*exp(-1i*beta);

output.r=r123;
output.t=t123;
output.r12=r12;
output.r23=r23;
output.t12=t12;
output.t23=t23;
output.t_single=t_single;
output.beta=beta;
