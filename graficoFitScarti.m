% GRAFICO FIT SCARTI
% Funzione che sfrutta linearFit per generare in automatico un grafico con
% pendenze intercetta e scarti. Comodo per un'analisi veloce dei set di
% dati raccolti.

function [A, B, dA, dB] = graficoFitScarti(data_x, data_y, sigma_x, sigma_y,fileName)
    arguments
        data_x (1, :) % vettore riga
        data_y (1, :) % vettore riga
        sigma_x (1, :) % vettore riga
        sigma_y (1, :) % vettore riga
        fileName (1, :) % vettore riga
    end
    
    % Variabili di setup ------------------
    
    % Unità di misura mostrate in legenda
    %unit_a ="F";
    %unit_pend = "F/m";

    unit_a ="F";
    unit_pend = "F/m";

    %y_label = "Capacità [F]";
    %x_label = "Lunghezza filo [m]";

    y_label = "Capacità [F]";
    x_label = "Lunghezza filo [m]";

    %t_itle = "Induttanza-Lunghezza";
    t_itle = "Capacità-Lunghezza";

    %-------------------------------------
    
    % Ricava a, b e sigma_totale da rappresentare sul grafico scarti.
    % Valori ottenuti con relative incertezze e chi.
    [a ,b, sa, sb, chi2fit] = linearFit(data_x, data_y, sigma_x, sigma_y);
    A=a; B=b; dA=sa; dB=sb;
    sigma_tot_vec = sqrt(sigma_y.^2 + (b*sigma_x).^2);
    %disp(sigma_tot_vec);
    
    % Impostazione layout
    figure();
    axes();
    h = 3;
    tl = tiledlayout(h,1);
    tl.TileSpacing = "tight";
    nexttile(1, [h-1, 1]);
    box on

   
    % Disegna retta regressione
    delta_x = max(data_x) - min(data_x);
    xlim([min(data_x)-0.1*delta_x max(data_x)+0.1*delta_x]);
    x2 = [min(data_x)-0.1*delta_x max(data_x)+0.1*delta_x];
    y2 = b*x2 + a;
    line(x2,y2,'Color','red','LineStyle','-')
    
    grid on;
    grid minor;
    hold on;

    % Plot dati x/y
    scatter(data_x, data_y, "MarkerEdgeColor",[0.00 0.45 0.74]);
    
    title(t_itle);
    ylabel(y_label);
    % textbox

    % Cifre significative da mostrare in legenda

    ta = numberToText(a, sa);
    tb = numberToText(b, sb);
    dof = (length(data_x)-2); 
    text = ["\alpha_L = " + ta + " " + unit_a; "\beta_L = " + tb + " " + unit_pend; "\chi_2 = " + fix(chi2fit) + "/" + dof];
    %text = ["\alpha = " + ta + " " + unit_a; "\beta = " + tb + " " + unit_pend];

    % Aggiungi legenda sopra o sotto la retta di fit in funzione del segno
    % della pendenza
    if (b < 0)
    annotation("textbox", [0.65,0.80,0.1,0.1], ...
        "BackgroundColor", [1,1,1], ...
        "FontSize", 14, ...
        "String", text ...
    );
    else
    annotation("textbox", [0.55,0.55,0.1,0.1], ...
        "BackgroundColor", [1,1,1], ...
        "FontSize", 14, ...
        "String", text ...
    );
    end  
    set(gca,'XTickLabel',[])
    set(gca, "FontSize", 14);
    
    % grafico scarti

    nexttile([1 1]);
    box on
    scarto_y = data_y - (b*data_x + a);
    e2 = errorbar(data_x, scarto_y, sigma_tot_vec);
    e2.LineStyle = 'none';
    xlim([min(data_x)-0.1*delta_x max(data_x)+0.1*delta_x]);
    ylim([-max(scarto_y.*2)-max(sigma_tot_vec)*1.7 max(scarto_y.*2)+max(sigma_tot_vec)*1.7]);
    x = [min(data_x)-0.1*delta_x max(data_x)+0.1*delta_x];
    y = [0 0];
    line(x,y,'Color','red','LineStyle','-')
    grid on;
    grid minor;
    hold on;
    scatter(data_x, scarto_y, "MarkerEdgeColor",[0.00 0.45 0.74]);
    ylabel("Scarto");
    xlabel(x_label);
    
    % Export figura in formato .png
    exportFigure(gcf, gca, fileName);
   
end