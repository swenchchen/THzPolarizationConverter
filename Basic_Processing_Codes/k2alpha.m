function alpha = k2alpha(k, freq)

c = 2.99792458; % speed of light 
alpha = 400*pi*freq.*k/c;  % in unit of cm-1