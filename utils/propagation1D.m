% -------------------------------------------------
% Funzione per propagare le incertezze su un modello
% a singolo parametro. 
% ---------------------------------------------------

% TODO
% Rendere l'intervallo infinitesimo dipendente dalla x della derivata

function sigma = propagation1D(model, data, sdata, h)
    arguments
        model (1,1) function_handle
        data (:,1) double {mustBeReal, mustBeFinite}
        sdata (:,1) double {mustBeReal, mustBeFinite}
        h (:,1) double {mustBeReal, mustBeFinite, mustBePositive} = 1e8;
    end

    % Derivata numerica del modello per propagazione incertezze
    d = (max(data) - min(data)) / h;
    d_model = @(x, dx) (model(x + dx) - model(x))/dx; 
    slope = d_model(data, d);
    sigma = slope.*sdata;
end