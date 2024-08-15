function output = unwrap_TDS (freq, fd_data, freq_range, phase_raw)

if nargin < 4
    phase_raw = 1;
end

fit_range = find(freq >freq_range(1) & freq<freq_range(2));
phase_unwrap_raw=unwrap(angle(fd_data));
poly_coef = polyfit(freq(fit_range), phase_unwrap_raw(fit_range), 1);
poly_shift = poly_coef(2); %take absolute shift
num_pi = round(poly_shift ./ (2.*pi)); 
phase_shift = num_pi .* 2 .* pi;
phase_unwrap=phase_unwrap_raw-phase_shift;

if phase_raw ==1
    output = phase_unwrap;
elseif phase_raw ==2 % output the phase shift at zero frequency
    output = phase_unwrap-poly_coef(1)*freq;
else % outpu the phase gradient
    output = poly_coef(1);
end