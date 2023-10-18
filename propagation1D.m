% -------------------------------------------------
% Funzione per propagare le incertezze su un modello
% a singolo parametro. 
% ---------------------------------------------------

function sigma = propagation1D(model, data, sdata)
    arguments
        model (1,1) function_handle
        data (:,1) double {mustBeReal, mustBeFinite}
        sdata (:,1) double {mustBeReal, mustBeFinite}
    end

    % Derivata numerica del modello per propagazione incertezze
    d = (max(data) - min(data)) / 1e8;
    d_model = @(x, dx) (model(x + dx) - model(x))/dx; 
    slope = d_model(data, d);
    sigma = slope.*sdata;
end