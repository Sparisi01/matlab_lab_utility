function text = numberToText(x, sx, cifre)
    % -------------------------------------------------
    % Funzione per la conversione di (x pm sx)
    % in testo. Le cifre significative  di x vengono 
    % aggiustate in automatico basandosi sul numero di 
    % cifre di sx.
    % x è un numero
    % sx è la sua incertezza
    % cifre (opzionale) è il numero di cifre significative 
    % -------------------------------------------------
    

    og = floor(log10(abs(x))); % ordine di grandezza
    % mantissa_x = round(x / (10^og) * 10^cifre) / (10^cifre);

    ogs = floor(log10(abs(sx)));
    sx = round(sx / (10^ogs)) * (10^ogs);

    if (~exist('cifre', 'var')) || cifre <= 0
        % correttamente uso la versione arrotondata di sx
        cifre = max(1, 1 + og - floor(log10(abs(sx))));
    end

    cifre = cifre - 1;

    mantissa_x = sprintf("%0." + cifre +  "f", x / (10^og));
    mantissa_sx = sx / (10^og) + "";

    if og == 0
        text = "(" + mantissa_x + " \pm " + mantissa_sx + ")";
    else
        text = "(" + mantissa_x + " \pm " + mantissa_sx + ")\times{}10^{" + og +  "}";
    end
end