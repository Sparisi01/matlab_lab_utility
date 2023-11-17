function [y_ricostruzione,s_y_ricostruzione] = ricostruisci_interpolazione(x_campioni, y_campioni, t , s_y_campioni)
    n_campioni = length(x_campioni);
    y_ricostruzione = zeros(size(t));
    s_y_ricostruzione = zeros(size(t));

    

    for ii = 1:n_campioni-1

        T = x_campioni(ii+1)-x_campioni(ii);

        y_ricostruzione = y_ricostruzione + ...
            (y_campioni(ii) * ((T - abs(t - x_campioni(ii)))/T .* (abs(t - x_campioni(ii)) < T)));   

        s_y_ricostruzione = s_y_ricostruzione + (...
            (s_y_campioni(ii) * ((T - abs(t - x_campioni(ii)))/T .* (abs(t - x_campioni(ii)) < T)))).^2;   
    end
    
    s_y_ricostruzione = sqrt(s_y_ricostruzione);
end