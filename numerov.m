% ---------------------------------------------------
% Funzione per la risoluzione numerica di equazioni 
% differenziali del secondo ordine nella forma:
% y(x)'' * (k(x)^2)*y(x) - s(x) = 0
% ---------------------------------------------------

function [X,Y] = numerov(Y_0, V_0, limits, n_steps, k, s)
    arguments
        Y_0 (1,1) double
        V_0 (1,1) double
        limits (1,2) double
        n_steps (1,1) 
        k (1,1) function_handle
        s (1,1) function_handle
    end

    % Limiti asse x
    x_min = limits(1);
    x_max = limits(2);

    if (x_min >= x_max)
        error("Limits must be two increasing values");
    end

    % Risoluzione asse x
    step = (x_max - x_min) / n_steps;

    % Vettore coordinate
    X = (x_min:step:x_max);
    Y = ones(size(X));

    % Condizioni iniziali
    Y(1) = Y_0;
    Y(2) = step*V_0 + Y_0;

    % Numerov method
    for j = 2:n_steps
        Y(j + 1) = ((s(X(j+1))+s(X(j-1))+10*s(X(j)))*step^2/12 + Y(j) * (2 -5/6 * step ^ 2 * k(X(j) ^ 2)) - Y(j - 1) * (1 + step^2/12 * k(X(j - 1)) ^ 2)) / (1 + step ^ 2/12 * k(X(j + 1)) ^ 2);
    end

end