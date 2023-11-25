function [y_rec,s_y_rec] = sincReconstruction(x_samples, y_samples, t , s_y_samples)
    n_samples = length(x_samples);
    y_rec = zeros(size(t));
    s_y_rec = zeros(size(t));

    for ii = 1:n_samples-1
        y_rec = y_rec + ...
            (y_samples(ii)*sinc(pi*(t - x_samples(ii))/((x_samples(ii+1)-x_samples(ii))*pi)));   

        s_y_rec = s_y_rec + (...
            (s_y_samples(ii)*sinc(pi*(t - x_samples(ii))/((x_samples(ii+1)-x_samples(ii))*pi)))).^2;   
    end
    
    s_y_rec = sqrt(s_y_rec);
end