function sigma = propagation(model, data, sdata)
    % Derivata numerica del modello per propagazione incertezze
    d = (max(data) - min(data)) / 1e8;
    slope = (model(data + d) - model(data))/d;
    sigma = slope.*sdata;
end