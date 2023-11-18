
close all;

t = 0:0.001:10;
y = 5*sin(t*2*pi*100) + 1*sin(t*2*pi*90);


ff = fourierTransform();
ff.data = y;
ff.dt = 0.001;



[frequencies, phases, amps,sigmaAmps,sigmaPhases, fig, ax] = ff.plotAbsTransform();

[peak_mean, peak_sigma, peak_amp, peak_amp_sigma, peak_phase, peak_phase_sigma] = ff.peakDetection();


fig;
hold on;
xlim(ax, [80,120]);
stem(peak_mean,peak_amp);

set(fig, 'visible', 'on'); 


