function [x_campioni, y_campioni] = campiona(t,V, f_campionamento, sfasamento)
    dt = t(2) - t(1);
    t_campionameto = 1/f_campionamento;
    n_campioni = floor(dt*length(t)/t_campionameto);
    x_campioni = (1:n_campioni)*t_campionameto + t(1) + sfasamento;
    y_campioni_index = floor((1:n_campioni)*t_campionameto/dt);
    y_campioni = V(y_campioni_index(:));
end