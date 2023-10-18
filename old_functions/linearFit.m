% REGRESSIONE LINEARE
% Funzione che esegue regressione linare con metodo iterativo.

function [res_a, res_b, err_a, err_b, chi] = linearFit(data_x, data_y, sigma_x, sigma_y)
    arguments
        data_x (1, :) % vettore riga
        data_y (1, :) % vettore riga
        sigma_x (1, :) % vettore riga
        sigma_y (1, :) % vettore riga
    end

    if ( ...
        length(data_x) ~= length(data_y) || ...
        length(data_x) ~= length(sigma_y) || ...
        length(data_x) ~= length(sigma_x) ...
    ) 
        throw("data_x, data_y, sigma_x and sigma_y should have the same length");
    end
    dataN = length(data_x);

    old_b = 9999999;
    a = 0;
    b = 0;
    togo = true;
    iterations = 0;
    max_it = 20;
    while(togo && iterations <= max_it)
        iterations = iterations + 1;
        new_old_b = b;

        % 1) usiamo b per calcolare un sigma_tot comprensivo di sigma_x e sigma_y
        sigma_tot = sqrt(sigma_y.^2 + (b*sigma_x).^2);
        sigma_tot2 = sigma_tot.^2;

        % 2) utilizziamo sigma_tot per calcolare a, b, sigma_a, sigma_b
        delta = (...
            sum(1./sigma_tot2)*sum(data_x.^2 ./ (sigma_tot2)) - ...
            (sum(data_x ./ (sigma_tot2)))^2 ...
        );
        a = (1 / delta) * ( ...
                sum(data_x.^2 ./ sigma_tot2)*sum(data_y ./ sigma_tot2) - ...
                sum(data_x ./ sigma_tot2)*sum(data_x.*data_y ./ sigma_tot2)...
            );
        b = (1 / delta) * ( ...
            sum(1./sigma_tot2)*sum(data_x.*data_y ./ sigma_tot2) - ...
            sum(data_x ./ sigma_tot2)*sum(data_y ./ sigma_tot2)...
        );
        sigma_a = sqrt(sum(data_x.^2 ./ sigma_tot2)/delta);
        sigma_b = sqrt(sum(1./sigma_tot2)/delta);

        delta_b = abs(old_b - b);
        old_b = new_old_b;

        %disp(" ");
        %disp("a: " + a);
        %disp("b: " + b);
        %disp("sigma_a: " + sigma_a);
        %disp("sigma_b: " + sigma_b);

        % 3) se delta_b è circa sigma_b ho finito altrimenti
        % disp("delta_b: " + delta_b);
        if(delta_b < sigma_b)
            togo = false;
        end

        % 1) ripetiamo l'operazione con il nuovo b
    end

    % Test del chi quadro
    sigma_tot = sqrt(sigma_y.^2 + (b*sigma_x).^2);
    chi_quadro = sum((data_y - data_x*b - a).^2 ./ (sigma_tot.^2));
    
    disp("-------- y = bx + a --------");
    disp("Coefficiente angolare (b): " + b + " ± " + sigma_b);
    disp("Reciproco (1/b): " + 1/b);
    disp("Intercetta (a): " + a + " ± " + sigma_a);
    disp("Chi2: " + chi_quadro + "/" + (length(data_x) - 2));
    disp("Iterazioni: " + iterations);
    
    expected_sigma_y = 1/(length(data_x) - 2) * sum((data_y - data_x*b - a).^2);
    disp("SigmaY attesa: " + expected_sigma_y);
    disp("SigmaY media: " + mean(sigma_y));
    disp("SigmaX media: " + mean(sigma_x));
    
    err_a = sigma_a;
    err_b = sigma_b;
    res_a = a;
    res_b = b;
    chi = chi_quadro;

end