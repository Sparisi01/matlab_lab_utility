function [x_samples, y_samples] = sampling(t,sampling_signal,f_sampling, offset)
    dt = t(2) - t(1);
    t_sampling = 1/f_sampling;
    n_samples = floor(dt*length(t)/t_sampling);
    x_samples = (1:n_samples)*t_sampling + t(1) + offset;
    y_samples_index = floor((1:n_samples)*t_sampling/dt);
    y_samples = sampling_signal(y_samples_index(:));
end