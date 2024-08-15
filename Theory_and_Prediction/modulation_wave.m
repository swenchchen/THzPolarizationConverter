t=0:0.1:140;
odd_ind = rem(round(t*2),2);

sq_wave = -30*ones(1,numel(t));
sq_wave(odd_ind==1) = 30;

mod_wave = ones(1, numel(t));
mod_wave(t>=0 & t<10) = 0;
mod_wave(t>=30 & t<110) = 0;
mod_wave(t>=130) = 0;

sq_wave_mod = mod_wave.*sq_wave;


figure
subplot(3,1,1)
plot(t, sq_wave)
ylim([-40,40]);
xlabel('time (ms)'); ylabel('Voltage (V)')
set(gca,'fontsize', 20)
subplot(3,1,2)
plot(t, 5*mod_wave)
ylim([-1,6]);
xlabel('time (ms)'); ylabel('TTL (V)')
set(gca,'fontsize', 20)
subplot(3,1,3)
plot(t, sq_wave_mod)
ylim([-40,40]);
xlabel('time (ms)'); ylabel('Voltage (V)')
set(gca,'fontsize', 20)