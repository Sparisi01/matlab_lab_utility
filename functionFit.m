% ---------------------------------------------------
% Funzione per eseguire fit a un modello generico e
% produrre grafici con numero di parametri variabile. 
% ---------------------------------------------------
% DIPENDENZE:
% - ./utils/numberToText.m
% - ./utils/exportFigure.m
% - ./utils/propagation1D.m
% ---------------------------------------------------

% TODO
% Finire implementazione showDataArray

classdef functionFit < handle

    properties
        datax (:, 1) double {mustBeReal, mustBeFinite}
        datay (:, 1) double {mustBeReal, mustBeFinite}
        sigmax (:, 1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        sigmay (:, 1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
        model (1, 1) function_handle = @(par, x) par(1) + x * par(2)
        par (:, 1) double {mustBeFinite}
        errpar (:, 1) double
        previousPar (:, 1) double {mustBeFinite}
        upperBounds (:, 1) double
        lowerBounds (:, 1) double
        units (:, 1) string
        parnames (:, 1) string
        yfit (:, 1) double {mustBeReal, mustBeFinite}
        chi2norm (1, 1) double {mustBeNonnegative}
        fig (1, 1)
        axes (:, 1)
        dof (1, 1) double {mustBeNonnegative}
        pValue (1, 1) double {mustBeNonnegative}
        showParArray (:, 1) logical
        showDataArray (:, 1) logical
        showChi (1, 1) logical
        showChiNorm (1, 1) logical
        showPValue (1, 1) logical
        showBox (1, 1) logical
        showScarti (1, 1) logical
        showGrid (1, 1) logical
        showZoom (1, 1) logical
        showInitialParModel (1,1) logical 
        showModel (1,1) logical
        continuosData (1,1) logical
        zoomPosition (1, 4) double {mustBeReal, mustBeFinite}
        modelColor (1, 4) double {mustBeReal, mustBeFinite}
        modelLineStyle (1, 1) string
        dataColor (1, 3) double {mustBeReal, mustBeFinite}
        lineWidth (1, 1) double
        name (1, 1) string
        labelx (1, 1) string
        labely (1, 1) string
        reslabely (1, 1) string
        logX (1, 1) logical
        logY (1, 1) logical
        xlim (1, 2) double {mustBeReal, mustBeFinite}
        ylim (1, 2) double {mustBeReal, mustBeFinite}
        resylim (1, 2) double {mustBeReal, mustBeFinite}
        boxPosition (1, 2) double {mustBeReal, mustBeFinite}
        pedice (1, 1) char
        fontSize (1, 1) double {mustBeReal, mustBeFinite}
        figureWidth (1, 1) double {mustBeReal, mustBeFinite}
        figureHeight (1, 1) double {mustBeReal, mustBeFinite}
        verbose (1, :) logical
        nDifferentialSteps (1, 1) double
    end

    methods

        % Valori di default per le opzioni
        function this = functionFit()
            % Dati --------------------------
            this.datax = [];
            this.datay = [];
            % Se le incertezze non vengono definite vengono inizializzate a
            % 1% dei dati. Se viene passato uno scalare la stessa
            % incertezza viene applicata a ogni punto.
            this.sigmax = [];
            this.sigmay = [];

            % Parametri modello -------------
            this.model = @(par, x) par(1) + x * par(2); % Modello su cui eseguire il fit
            this.par = []; % Valore dei parametri
            this.previousPar = []; % Valore dei parametri
            this.errpar = []; % Errore dei parametri
            this.upperBounds = []; % UpperBound parametri
            this.lowerBounds = []; % LowerBound parametri
            this.units = []; % Units for parameters in legend box
            this.parnames = []; % Parameters name in legend box
            
            % Risultati fit  ----------------
            this.yfit = []; % Y calcolati post regressione con parametri ottimizzati
            this.chi2norm = inf; % CHi quadro normalizzato fit
            this.dof = 0; % Gradi di libertà fit
            this.pValue = 0; % P value fit
            this.fig = 0; % Fig dopo aver generato figura e residui
            this.axes = []; % Array contenenti gli assi dopo aver generato figura e residui
            
            % Estetica ----------------------
            this.showParArray = [1 1]; % Scegli quali parametri mostrate con un array di bool della stessa dimensione di par
            this.showDataArray = []; % Scegli quali dati mostrate con un array di bool della stessa dimensione dei dati. I dati non mostrati in grafico contribuiscono comunque al fit
            this.showChi = 1; % Mostra o no chi quadro in legenda
            this.showChiNorm = 0; % Mostra chi2normalizzato
            this.showPValue = 0; % Mostra PValue
            this.showBox = 1; % Mostra o no box parametri
            this.showScarti = 1; % Mostra o no grafico degli scarti
            this.showInitialParModel = 0; % Mostra modello con parametri iniziali su grafico
            this.showModel = 1; % Mostra modello sul grafico
            this.continuosData = 0; % Modalità visualizzazione continua
            this.modelColor = [1 0 0 1]; % Colore linea modello
            this.modelLineStyle = '-'; % Stile linea modello e retta scarti
            %this.dataColor = [0.00 0.45 0.74]; % Colore dati bello
            this.dataColor = [0.00 0.00 1.00]; % Colore dati blu
            this.lineWidth = 2; % Spessore linea modello
            this.xlim = [0 0]; % Xlim, se uguali o in ordine sbagliato viene impostato in automatico
            this.ylim = [0 0]; % Ylim, se uguali o in ordine sbagliato viene impostato in automatico
            this.resylim = [0 0]; % Ylim residui, se uguali o in ordine sbagliato viene impostato in automatico
            this.verbose = 1; % Bla Bla Bla
            this.name = "Model"; % Titolo grafico
            this.labelx = "X Axes"; % Label Asse x
            this.labely = "Y Axes"; % Label Asse y
            this.reslabely = "Scarti"; % Label Asse y scarti
            this.logX = 0; % Asse X logaritmico (anche per scarti garantendo allineamento)
            this.logY = 0; % Asse Y logaritmico
            this.boxPosition = [0.55, 0.55]; % [x, y] la dimensione di aggiusta in automatico
            this.pedice = ' '; % Pedice parametri legenda. Utile se si hanno molti grafici con parametri omonomi.
            this.showZoom = false; % Mostra grafico con zoom su un punto e barre incertezza
            this.zoomPosition = [0.21, 0.75, 0.15, 0.15]; % [x, y, w, h]
            this.showGrid = 1; % Mostra griglia minor sui grafici
            this.fontSize = 14; % Dimensione font sia nelle label che nella box
            this.figureWidth = 8; % Larghezza immagine salvata in pollici
            this.figureHeight = 6; % Altezza immagine salvata in pollici            
        end

        % Genera immagine plot usando il modello dato
        function [par, errpar, yfit, chi2norm, dof, pValue, fig, ax] = plotModelFit(this, file_name, showFig)
            arguments
                this
                file_name (1, 1) string = "",
                showFig (1,1) logical = 1
            end
            
            % Salva parametri precedenti al fit
            this.previousPar = this.par;

            % Esegui fit lineare e salva negli attributi dell'oggetto i nuovi valori dei parametri
            [par, errpar, yfit, chi2norm, dof, pValue] = modelFit(this);
                       
            % Genera Figura
            [fig, ax] = generatePlotFig(this, false);
            this.fig = fig;
            this.axes = ax;

            % Rendi visibile figura
            if showFig
                set(fig, 'visible', 'on'); 
            end

            % Esporta figura in formato png
            if (strlength(file_name) > 0)
                exportFigure(fig, ax, file_name, this.fontSize, this.figureWidth, this.figureHeight);
            end

        end

        % Genera immagine plot usando il modello lineare
        function [par, errpar, yfit, chi2norm, dof, pValue, fig, ax] = plotLinearFit(this, file_name, showFig)
            arguments
                this
                file_name (1, 1) string = "",
                showFig (1,1) logical = 1
            end
            
            % Salva parametri precedenti al fit
            this.previousPar = this.par;

            % Esegui fit non lineare e salva negli attributi dell'oggetto i nuovi valori dei parametri
            [par, errpar, yfit, chi2norm, dof, pValue] = linearFit(this);
                                   
            % Genera Figura
            [fig, ax] = generatePlotFig(this, true);
            this.fig = fig;
            this.axes = ax;
            
            % Rendi visibile figura
            if showFig
                set(fig, 'visible', 'on'); 
            end
            
            % Esporta figura in formato png
            if (strlength(file_name) > 0)
                exportFigure(fig, ax, file_name, this.fontSize, this.figureWidth, this.figureHeight);
            end

        end

        % Funzione per fit lineare
        function [par, errpar, yfit, chi2norm, dof, pValue] = linearFit(this)

            % Contolla validità parametri
            safetyCheck(this)

            data_x = this.datax;
            data_y = this.datay;
            sigma_x = this.sigmax;
            sigma_y = this.sigmay;
            
            dataN = length(data_x);
            old_b = 9999999;
            a = 0;
            b = 0;
            togo = true;
            iterations = 0;
            max_it = 20;

            while (togo && iterations <= max_it)
                iterations = iterations + 1;
                new_old_b = b;

                % 1) usiamo b per calcolare un sigma_tot comprensivo di sigma_x e sigma_y
                sigma_tot = sqrt(sigma_y .^ 2 + (b * sigma_x) .^ 2);
                sigma_tot2 = sigma_tot .^ 2;

                % 2) utilizziamo sigma_tot per calcolare a, b, sigma_a, sigma_b
                delta = ( ...
                    sum(1 ./ sigma_tot2) * sum(data_x .^ 2 ./ (sigma_tot2)) - ...
                    (sum(data_x ./ (sigma_tot2))) ^ 2 ...
                );
                a = (1 / delta) * ( ...
                    sum(data_x .^ 2 ./ sigma_tot2) * sum(data_y ./ sigma_tot2) - ...
                    sum(data_x ./ sigma_tot2) * sum(data_x .* data_y ./ sigma_tot2) ...
                );
                b = (1 / delta) * ( ...
                    sum(1 ./ sigma_tot2) * sum(data_x .* data_y ./ sigma_tot2) - ...
                    sum(data_x ./ sigma_tot2) * sum(data_y ./ sigma_tot2) ...
                );
                sigma_a = sqrt(sum(data_x .^ 2 ./ sigma_tot2) / delta);
                sigma_b = sqrt(sum(1 ./ sigma_tot2) / delta);

                delta_b = abs(old_b - b);
                old_b = new_old_b;

                % 3) ripetiamo l'operazione con il nuovo b se la variazione
                % Ã¨ minore dell'incertezza sul b
                if (delta_b < sigma_b)
                    togo = false;
                end

            end

            % Test del chi quadro
            sigma_tot = sqrt(sigma_y .^ 2 + (b * sigma_x) .^ 2);
            chi2 = sum((data_y - data_x * b - a) .^ 2 ./ (sigma_tot .^ 2));
            
            % Results
            yfit = data_x * b + a;
            par(1) = a;
            par(2) = b;
            errpar(1) = sigma_a;
            errpar(2) = sigma_b;
            chi2norm = chi2 / (dataN - 2);
            dof = dataN - 2;
            pValue = chi2cdf(chi2, dof);
            
            this.par = par;
            this.errpar = errpar;
            this.chi2norm = chi2norm;
            this.yfit = yfit;
            this.dof = dof;
            this.pValue = pValue;

        end

        % Funzione per fit non lineare
        function [par, errpar, yfit, chi2norm, dof, pValue, flag] = modelFit(this)

            % Contolla validità parametri
            safetyCheck(this)

            % Definizione funzione scarti a partire dal modello
            scarti = @(par, xd, yd, ed) (this.model(par, xd) - yd) ./ ed;

            % Aggiusta dimensioni upper e lower bound. Se non impostate
            % vengono inizializzati a due vettori delle dimensione di par a
            % valori + e - inf
            if (isempty(this.upperBounds))
                tmp_ub = ones(size(this.par)) * inf;
            else
                tmp_ub = this.upperBounds;
            end

            if (isempty(this.lowerBounds))
                tmp_lb = ones(size(this.par)) * -inf;
            else
                tmp_lb = this.lowerBounds;
            end
           
            if this.verbose
                options = optimoptions('lsqnonlin');
            else
                 options = optimoptions('lsqnonlin', 'Display', 'none');
            end
            
            % Fit non lineare Libreria Matlab
            [par, resnorm, ~, flag, ~, ~, jacobian] = lsqnonlin(scarti, this.par, tmp_lb, tmp_ub, options, this.datax, this.datay, this.sigmay);
                                 
            % Gradi di liberta
            dof = (length(this.datax) - length(par));
            
            % Yfit ottimizzati
            yfit = this.model(par, this.datax);
            
            % Propaga incertezze x lungo y
            sigma_propagata = propagation1D(@(x) this.model(par, x), this.datax, this.sigmax);
            sigmaScarti = sqrt(this.sigmay .^2 + (sigma_propagata).^2); 
            
            % Calcolo Chi quadro e pValue
            scarto_y = this.datay - yfit;            
            chi2norm = sum((scarto_y./sigmaScarti).^2)/dof;                      
            pValue = 1 - chi2cdf(resnorm, dof);    

            % Matrice di covariaza
            covar = inv(jacobian' * jacobian);

            % Calcolo errori parametri
            var = diag(covar);
            sigma = sqrt(var);
            sigmaf = full(sigma);
            errpar = sigmaf * sqrt(chi2norm);           
                        
            % Assegna variabili oggetto
            this.par = par;
            this.errpar = errpar;
            this.chi2norm = chi2norm;
            this.yfit = yfit;
            this.dof = dof;
            this.pValue = pValue;            
        end
    end

    methods (Hidden)

        % Fa tante cose belle
        function safetyCheck(this)
            
            % Controlla dimensione dati -----------------------------------
            if any(size(this.datay) ~= size(this.datax))                   
                error('u:stuffed:it', 'datax, datay devono essere della stesa dimensione');
            end

            % Controlla dimensione bound ----------------------------------
            if ~(isempty(this.upperBounds) || isempty(this.lowerBounds)) && ...
                    any(size(this.upperBounds) ~= size(this.lowerBounds))

                error('u:stuffed:it', "ub e lb devono essere della stesa dimensione");
            end
            
            % Controlla dimensione incertezze -----------------------------
            if (isempty(this.sigmax))
                this.sigmax = zeros(size(this.datay));
                warning("Incertezze su datax non assegnate. Inizializzate a 0")           
            else
                if size(this.sigmax) == [1,1]
                    this.sigmax = ones(size(this.datax)) * this.sigmax;
                end                
            end
            
            % Controlla dimensioni sigma e inizializza ---------------------
            if (isempty(this.sigmay))
                this.sigmay = ones(size(this.datay));
                warning("Incertezze su datay non assegnate. Inizializzate a un vettore unitario.")        
            else
                if size(this.sigmay) == [1,1]
                    this.sigmay = ones(size(this.datay)) * this.sigmay;
                end 
            end
        end

        function [fig, ax] = generatePlotFig(this, isLinearFit)

            fig = figure('units','inch','position',[0,0,this.figureWidth,this.figureHeight],'visible','off');
            
            % Costrusci tiled layout
            if this.showScarti
                h = 3;
                tiledlayout(h, 1);               
                ax(1) = nexttile(1, [h - 1, 1]);
            else
                ax(1) = axes();
            end
            
            % Determina se i dati vanno filtrati
            toFilter = 1;
            if (isempty(this.showDataArray))
                this.showDataArray = ones(size(this.datax));
                toFilter = 0;
            else
                if size(this.showDataArray) ~= size(this.datax)
                    error('u:stuffed:it', "showDataArray e datax, datay devono essere della stessa dimensione.");
                end
            end
                    
            % Filtra dati in base a showDataArray
            if toFilter
                jj = 1;
                filtered_datax = ones(sum(this.showDataArray),1);
                filtered_datay = ones(sum(this.showDataArray),1);

                for ii = 1:length(this.showDataArray)
                    if this.showDataArray(ii) == 1                        
                        filtered_datax(jj) = this.datax(jj);
                        filtered_datay(jj) = this.datay(jj);
                        jj = jj + 1;
                    end
                end
            else
                filtered_datax = this.datax;
                filtered_datay = this.datay;          
            end            
                        
            hold on
            box on

            % Plot dei dati
            scatter(filtered_datax, filtered_datay, "MarkerEdgeColor", this.dataColor);
                        
            % Disegna la funzione con parametri ottimizzati e parametri iniziali     
            if isLinearFit
                if this.showModel
                    fplot(@(x) this.par(1) + x * this.par(2), 'Color', this.modelColor, 'LineStyle', this.modelLineStyle, "LineWidth",this.lineWidth);
                end
                if this.showInitialParModel
                    hold on;
                    fplot(@(x) this.previousPar(1) + x * this.previousPar(2), 'Color', 'k', 'LineStyle', this.modelLineStyle,"LineWidth",this.lineWidth);
                end             
            else
                if this.showModel
                    fplot(@(x) this.model(this.par, x), 'Color', this.modelColor, 'LineStyle', this.modelLineStyle,"LineWidth",this.lineWidth);
                end
                if this.showInitialParModel
                    hold on;
                    fplot(@(x) this.model(this.previousPar, x), 'Color', 'k', 'LineStyle', this.modelLineStyle,"LineWidth",this.lineWidth);
                end
            end

            % Imposta dinamicamente i limiti ------------------------------

            delta_x = abs(max(this.datax) - min(this.datax));
            delta_y = abs(max(this.datay) - min(this.datay));
            
            % Limiti asse X scarti
            if this.xlim(1) >= this.xlim(2)
                % Gestione caso asse logaritmico per evitare il passaggio
                % dei limiti attraverso lo 0
                if this.logX
                    xlim(sort([min(this.datax) * 0.2 max(this.datax) * 1.2]));
                else     
                    xlim(sort([min(this.datax) - 0.1 * delta_x max(this.datax) + 0.1 * delta_x]));
                end
            else
                xlim([this.xlim(1) this.xlim(2)]);
            end

            if this.ylim(1) >= this.ylim(2)  
                % Gestione caso asse logaritmico per evitare il passaggio
                % dei limiti attraverso lo 0
                if this.logY
                    ylim(sort([min(this.datay) * 0.2 max(this.datay) * 1.2]));
                else     
                    ylim(sort([min(this.datay) - 0.1 * delta_y max(this.datay) + 0.1 * delta_y]));
                end             
            else
                ylim([this.ylim(1) this.ylim(2)]);
            end
            
            % Estetica assi -----------------------------------------------

            title(this.name);
            ylabel(this.labely);
            
            if this.showGrid
                grid minor;
            end 

            if this.logX
                set(ax(1), 'XScale', 'log')
            end

            if this.logY
                set(ax(1), 'YScale', 'log')
            end

            if ~this.showScarti
                xlabel(this.labelx);
            else
                set(gca, 'XTickLabel', []);
            end
            
            % Box parametri -----------------------------------------------
            
            if this.showBox
                              
                % Aggiusta array unita di misura.
                if (isempty(this.units))       
                    for ii = 1:length(this.par)
                        this.units(ii) = "";
                    end                  
                else
                    if any(size(this.units) ~= size(this.par))
                        for ii = 1:length(this.par)
                            this.units(ii) = "";
                        end
                        warning("Dimensione di par diversa da units. Parametro inizializzato in automatico.")
                    end
                end
                
                % Aggiusta array nome parametri. Se non impostato viene inizializzato a un vettore delle
                % dimensioni di par formato da caratteri crescenti in ordine
                % alfabetico
                if (isempty(this.parnames))
                    for ii = 1:length(this.par)
                        this.parnames(ii) = 'a' + ii - 1;
                    end                
                else
                    if any(size(this.parnames) ~= size(this.par))
                        for ii = 1:length(this.par)
                            this.parnames(ii) = 'a' + ii - 1;
                        end
                        warning("Dimensione di par diversa da parnames. Parametro inizializzato in automatico.")
                    end
                end              
                
                if isempty(this.showParArray)
                    this.showParArray = ones(size(this.par));
                else
                    if any(size(this.showParArray) ~= size(this.par))
                        this.showParArray = ones(size(this.par));
                        warning("Dimensione di par diversa da showParArray. Parametro inizializzato in automatico.")
                    end
                end

                % Costruisci dinamicamente la legenda con n parametri.
                txt = "";            
                for ii = 1:length(this.par)
                    if this.showParArray(ii)
                        t = numberToText(this.par(ii), this.errpar(ii));                                            
                        if (this.pedice ~= ' ')
                            txt(length(txt) + 1) = this.parnames(ii) + "_{" + this.pedice + "} = " + t + " " + this.units(ii);
                        else
                            txt(length(txt) + 1) = this.parnames(ii) + " = " + t + " " + this.units(ii);
                        end
                    end
                end

                % Se il chi2 è minore di 2 vengono tenute 2 cifre significative
                % sennò si arrotonda all'intero.
                if this.chi2norm * this.dof > 2
                    nRound = 0;
                else
                    nRound = 2;
                end

                if this.showChi
                    txt(length(txt) + 1) = "\chi^2_{" + this.pedice + "} = " + round(this.chi2norm * this.dof, nRound) + "/" + this.dof;
                end

                if this.showChiNorm
                    txt(length(txt) + 1) = "\chi^2_{" + this.pedice + "} = " + round(this.chi2norm, 2);
                end

                if this.showPValue
                    txt(length(txt) + 1) = "pValue_{" + this.pedice + "} = " + round(this.pValue, 2);
                end

                % Posizione box
                annotation("textbox", [this.boxPosition 0 0], ...
                    "BackgroundColor", [1, 1, 1], ...
                    "fontSize", this.fontSize, ...
                    "String", txt(2:end), ...
                    'FitBoxToText', 'on' ...
                );
            end


            % Grafico scarti  ---------------------------------------------

            if this.showScarti                
                ax(2) = nexttile([1 1]);
               
                box on

                % Propagazione incertezze sugli scarti
                if isLinearFit
                    sigmaScarti = sqrt(this.sigmay .^ 2 + (this.par(2) * this.sigmax) .^ 2);
                else
                    % Derivata numerica del modello per propagazione incertezze                   
                    sigma_propagata = propagation1D(@(x) this.model(this.par, x), this.datax, this.sigmax);
                    sigmaScarti = sqrt(this.sigmay .^2 + (sigma_propagata).^2);    
                end
                
                % Vettore scarti
                scarto_y = this.datay - this.yfit;

                % Filtra scarti in base a showDataArray                           
                if toFilter
                    jj = 1;
                    filtered_scarto_y = ones(sum(this.showDataArray),1);
                    filtered_s_scarto_y = ones(sum(this.showDataArray),1);
                    for ii = 1:length(this.showDataArray)
                        if this.showDataArray(ii) == 1
                            filtered_scarto_y(jj) = scarto_y(jj);
                            filtered_s_scarto_y(jj) = sigmaScarti(jj);
                            jj = jj + 1;
                        end
                    end
                else
                    filtered_scarto_y = scarto_y;
                    filtered_s_scarto_y = sigmaScarti;
                end
                                
                hold on;
                                
                if this.showGrid
                    grid minor;
                end

                if this.logX
                    set(ax(2), 'XScale', 'log')
                end

                % Linea orizzontale lungo y=0 negli scarti
                % Gestione caso asse logaritmico per evitare il passaggio
                % delle coordinate attraverso lasse 0
                if this.logX                   
                    x = [sort([min(this.datax) * 0.2 max(this.datax) * 1.2])];
                    y = [0 0];
                else
                    x = [min(this.datax) - 0.1 * delta_x max(this.datax) + 0.1 * delta_x];
                    y = [0 0];
                end                
                line(x, y, 'Color', this.modelColor, 'LineStyle', '-',"LineWidth", this.lineWidth)
                
                % Plot degli scarti
                if this.continuosData
                    plot(filtered_datax, filtered_scarto_y,"LineWidth",this.lineWidth,"Color",[this.dataColor 1]);
                    plot(filtered_datax, filtered_scarto_y - filtered_s_scarto_y, "LineStyle","--","Color", [this.dataColor .3],"LineWidth",this.lineWidth);
                    plot(filtered_datax, filtered_scarto_y + filtered_s_scarto_y, "LineStyle","--","Color", [this.dataColor .3],"LineWidth",this.lineWidth);
                else
                    e2 = errorbar(filtered_datax, filtered_scarto_y, filtered_s_scarto_y);
                    e2.LineStyle = 'none';                  
                    scatter(this.datax.*this.showDataArray, scarto_y.*this.showDataArray, "MarkerEdgeColor", this.dataColor);
                end
               
                % Limiti asse X scarti
                if this.xlim(1) >= this.xlim(2)
                    % Gestione caso asse logaritmico per evitare passaggio
                    % dei limiti attraverso lo 0
                    if this.logX
                        xlim(sort([min(this.datax) * 0.2 max(this.datax) * 1.2]));
                    else     
                        xlim(sort([min(this.datax) - 0.1 * delta_x max(this.datax) + 0.1 * delta_x]));
                    end
                else
                    xlim([this.xlim(1) this.xlim(2)]);
                end
                
                % Limiti asse Y scarti
                if this.resylim(1) >= this.resylim(2)
                    tmp_limits = max([abs(max(scarto_y)) + 2*abs(max(sigmaScarti)) abs(min(scarto_y)) + 2*abs(max(sigmaScarti))]);
                    ylim([-tmp_limits tmp_limits]); % Asse Y scarti simmetrico rispetto allo 0      
                else
                    ylim([this.resylim(1) this.resylim(2)])
                end
                
                % Label Scarti
                ylabel(this.reslabely);
                xlabel(this.labelx);                                              
            end

            % Grafico errore zoom  ----------------------------------------

            if (this.showZoom)
                hold on;
                id = round(length(this.datax) / 2);
                x = this.datax(id);
                y = this.datay(id);
                ax(length(ax)+1) = axes("Position", this.zoomPosition);
                errorbar(x, y, -this.sigmay(id), this.sigmay(id), -this.sigmax(id), this.sigmax(id), 'o');
                xlim(ax(length(ax)), [(x - this.sigmax(id) * 1.5) (x + this.sigmax(id) * 1.5)]);
                ylim(ax(length(ax)), [(y - this.sigmay(id) * 1.5) (y + this.sigmay(id) * 1.5)]);
                if this.showGrid
                    grid on;
                end
            end

            for a = ax
                set(a, "FontSize", this.fontSize);
            end

        end        
    end
end
