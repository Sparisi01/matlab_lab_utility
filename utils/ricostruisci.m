function [y_ricostruzione] = ricostruisci(x_campioni, y_campioni, t)
    n_campioni = length(x_campioni);
    y_ricostruzione = zeros(size(t));
    for ii = 1:n_campioni-1
        y_ricostruzione = y_ricostruzione + ...
            (y_campioni(ii)*sinc(pi*(t - x_campioni(ii))/((x_campioni(ii+1)-x_campioni(ii))*pi)));   
    end
end