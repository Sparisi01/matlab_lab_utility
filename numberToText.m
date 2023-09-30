function text = numberToText(x, sx)
% x è un numero
% sx è la sua incertezza
% cifre è il numero di cifre significative

    og = floor(log10(abs(x))); % ordine di grandezza
    % mantissa_x = round(x / (10^og) * 10^cifre) / (10^cifre);

    ogs = floor(log10(abs(sx)));
    sx = round(sx / (10^ogs)) * (10^ogs);
    cifre = strlength(string(sx / (10^og)))-2;

    mantissa_x = sprintf("%0." + cifre +  "f", x / (10^og));
    mantissa_sx = sx / (10^og) + "";

    if og == 0
        text = "(" + mantissa_x + " \pm " + mantissa_sx + ")";
    else
        text = "(" + mantissa_x + " \pm " + mantissa_sx + ")\times{}10^{" + og +  "}";
    end
end