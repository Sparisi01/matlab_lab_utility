function [y_rec,s_y_rec] = triangularReconstruction(x_samples, y_samples, t , s_y_samples)
    n_samples = length(x_samples);
    y_rec = zeros(size(t));
    s_y_rec = zeros(size(t));

    

    for ii = 1:n_samples-1

        T = x_samples(ii+1)-x_samples(ii);

        y_rec = y_rec + ...
            (y_samples(ii) * ((T - abs(t - x_samples(ii)))/T .* (abs(t - x_samples(ii)) < T)));   

        s_y_rec = s_y_rec + (...
            (s_y_samples(ii) * ((T - abs(t - x_samples(ii)))/T .* (abs(t - x_samples(ii)) < T)))).^2;   
    end
    
    s_y_rec = sqrt(s_y_rec);
end