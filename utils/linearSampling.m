function ysamples = linearSampling(datax, datay, xsamples, varargin)
    % function ysamples = LinearSampling(datax, datay, xsamples)
    %
    % Fornite coppie di dati (datax, datay) restituisce i valori 'ysamples'
    % corrispondenti a entrate 'xsamples' interpolando linearmente tra i dati.
    % La funzione suppone che datax abbia valori in ordine crescente.


    % SAFETY CHECKS
    if length(datax) ~= length(datay)
        assert(false, "datax and datay should have the same length");
    end
    % if max(xsamples) > max(datax) | min(xsamples) < min(datay)
    %     disp("max(xsamples) = " + max(xsamples));
    %     disp("max(datax) = " + max(datax));
    %     disp("max(xsamples) > max(datax) = " + (max(xsamples) > max(datax)));
    %     disp("min(xsamples) = " + min(xsamples));
    %     disp("min(datax) = " + min(datax));
    %     disp("min(xsamples) < min(datay) = " + (min(xsamples) < min(datay)));
    %     error("xsamples should have values inside the interval [min(datax), max(datax)]");
    % end

    ascending = true;
    if nargin > 3 & varargin(4) == "descending"
        ascending = false;
    end

    % LINEAR INTERPOLATION (lerp)
    ysamples = [];
   
    for x = xsamples
        id_before = 1;
        id_after = 1;
        for kk = 1:length(datax)
            if ascending
                if(datax(kk) > x)
                    id_before = kk - 1;
                    id_after = kk;
                    break;
                end
            else
                if(datax(kk) < x)
                    id_before = kk - 1;
                    id_after = kk;
                    break;
                end
    
            end
        end
        y = 0;  
        if x >= max(datax)
            y = max(datax);
        elseif x <= min(datax)
            y = min(datax);
        else
            if(id_before == 0)
                id_before = id_before + 1;
                id_after = id_after + 1;
            end
            y = datay(id_before) + (datay(id_after) - datay(id_before))*(x - datax(id_before))/(datax(id_after) - datax(id_before));
        end
        ysamples = [ysamples, y];
    end
  end