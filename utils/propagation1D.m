% -------------------------------------------------
% Funzione per propagare le incertezze su un modello
% a singolo parametro. 
% ---------------------------------------------------

function sigma = propagation1D(model, data, sdata)
    arguments
        model (1,1) function_handle
        data (:,1) double {mustBeReal, mustBeFinite}
        sdata (:,1) double {mustBeReal, mustBeFinite, mustBeNonnegative}        
    end

    % Metodo differenziale variabile --------------------------------------
    sigma = zeros(size(data));
    d_sum = 0;

    % Differenziale modello
    d_model = @(x, dx) (model(x + dx) - model(x))/dx; 

    for ii = 1:length(data)
        if ii == length(data)
            d = d_sum/(ii-1);
        else
            d = data(ii+1) - data(ii);
            if d == 0
                if d_sum == 0
                    d = 1e-15;
                else
                    d = d_sum/ii;
                end
            end
            d_sum = d_sum + d;
        end       
        
        slope = d_model(data(ii), d);
        sigma(ii) = abs(slope.*sdata(ii));        
    end
end