classdef functionFit < handle
    % ---------------------------------------------------
    % Funzione per eseguire fit a un modello generico e
    % produrre grafici con numero di parametri variabile. 
    % ---------------------------------------------------
    % DIPENDENZE:
    % - ./utils/numberToText.m
    % - ./utils/exportFigure.m
    % - ./utils/propagation1D.m
    % - ./utils/avoidOversampling.m
    % - ./utils/tight_subplot.m
    % ---------------------------------------------------

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
        noOversampling (1,1) logical
        yfit (:, 1) double {mustBeReal, mustBeFinite}
        chi2norm (1, 1) double {mustBeNonnegative}
        fig (1, 1)
        axes (:, 1)
        dof (1, 1) double {mustBeNonnegative}
        pValue (1, 1) double {mustBeNonnegative}
        showParArray (:, 1) logical
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
        hMargins (1, 2) double {mustBeReal, mustBeFinite}
        wMargins (1, 2) double {mustBeReal, mustBeFinite}
        padding (1, 1) double {mustBeReal, mustBeFinite}
        verbose (1, :) logical
    end

    methods

        % Valori di default per le opzioni
        function this = functionFit()
            % Dati --------------------------
            this.datax = [];
            this.datay = [];
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
            this.noOversampling = 0;  % Merge data on the same x

            % Risultati fit  ----------------
            this.yfit = []; % Y calcolati post regressione con parametri ottimizzati
            this.chi2norm = inf; % CHi quadro normalizzato fit
            this.dof = 0; % Gradi di libertà fit
            this.pValue = 0; % P value fit
            this.fig = 0; % Fig dopo aver generato figura e residui
            this.axes = []; % Array contenenti gli assi dopo aver generato figura e residui
            
            % Estetica ----------------------
            this.showParArray = [1 1]; % Scegli quali parametri mostrate con un array di bool della stessa dimensione di par
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
            this.boxPosition = [0.11 0.89]; % [x, y] la dimensione di aggiusta in automatico
            this.pedice = ' '; % Pedice parametri legenda. Utile se si hanno molti grafici con parametri omonomi.
            this.showZoom = false; % Mostra grafico con zoom su un punto e barre incertezza
            this.zoomPosition = [0.70 0.70, 0.15, 0.15]; % [x, y, w, h]
            this.showGrid = 1; % Mostra griglia minor sui grafici
            this.fontSize = 14; % Dimensione font sia nelle label che nella box
            this.figureWidth = 8; % Larghezza immagine salvata in pollici
            this.figureHeight = 6; % Altezza immagine salvata in pollici    
            this.hMargins = [0.12 0.08];
            this.wMargins = [0.13 0.13]; 
            this.padding = 0.05; 
        end

        % Genera immagine plot usando il modello dato
        function [par, errpar, yfit, chi2norm, dof, pValue, fig, ax] = plotModelFit(this, file_name, showFig)
            arguments
                this
                file_name (1, 1) string = "",
                showFig (1,1) logical = 1
            end
            
            this.previousPar = this.par;

            [par, errpar, yfit, chi2norm, dof, pValue] = modelFit(this);
                       
            [fig, ax] = generatePlotFig(this, false);
            this.fig = fig;
            this.axes = ax;

            if showFig
                set(fig, 'visible', 'on'); 
            end

            if (strlength(file_name) > 0)
                exportFigure(fig, ax, file_name, this.fontSize, this.figureWidth, this.figureHeight);
            end

        end

        function [par, errpar, yfit, chi2norm, dof, pValue, fig, ax] = plotLinearFit(this, file_name, showFig)
            arguments
                this
                file_name (1, 1) string = "",
                showFig (1,1) logical = 1
            end
            
            this.previousPar = this.par;

            [par, errpar, yfit, chi2norm, dof, pValue] = linearFit(this);
                                   
            [fig, ax] = generatePlotFig(this, true);
            this.fig = fig;
            this.axes = ax;
            
            if showFig
                set(fig, 'visible', 'on'); 
            end
            
            if (strlength(file_name) > 0)
                exportFigure(fig, ax, file_name, this.fontSize, this.figureWidth, this.figureHeight);
            end

        end

        function [par, errpar, yfit, chi2norm, dof, pValue] = linearFit(this)

            safetyCheck(this);

            if(this.noOversampling)
                [this.datax, this.datay, this.sigmay] = avoidOversampling(this.datax, this.datay, this.sigmay);
            end

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

        function [par, errpar, yfit, chi2norm, dof, pValue, flag] = modelFit(this)
            
            % Contolla validità parametri
            safetyCheck(this)

            if(this.noOversampling)
                [this.datax, this.datay, this.sigmay] = avoidOversampling(this.datax, this.datay, this.sigmay);
            end           
            
            scarti = @(par, xd, yd, ed) (this.model(par, xd) - yd) ./ ed;

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
            
            [par, resnorm, ~, flag, ~, ~, jacobian] = lsqnonlin(scarti, this.par, tmp_lb, tmp_ub, options, this.datax, this.datay, this.sigmay);
                      
            dof = (length(this.datax) - length(par));  
            yfit = this.model(par, this.datax);
            
            sigma_propagata = propagation1D(@(x) this.model(par, x), this.datax, this.sigmax);
            sigmaScarti = sqrt(this.sigmay .^2 + (sigma_propagata).^2); 
            
            scarto_y = this.datay - yfit;            
            chi2norm = sum((scarto_y./sigmaScarti).^2)/dof;                      
            pValue = 1 - chi2cdf(resnorm, dof);    

            covar = inv(jacobian' * jacobian);
            var = diag(covar);
            sigma = sqrt(var);
            sigmaf = full(sigma);           
            errpar = sigmaf * sqrt(chi2norm);           
                        
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
                if length(this.sigmax) == 1
                    this.sigmax = ones(size(this.datax)) * this.sigmax;
                end                
            end
            
            % Controlla dimensioni sigma e inizializza ---------------------
            if (isempty(this.sigmay))
                this.sigmay = ones(size(this.datay));
                warning("Incertezze su datay non assegnate. Inizializzate a un vettore unitario.")        
            else
                if length(this.sigmay) == 1
                    this.sigmay = ones(size(this.datay)) * this.sigmay;
                end 
            end
        end

        function [fig,ax] = generatePlotFig(this, isLinearFit)

            fig = figure('units','inch','position',[0,0,this.figureWidth,this.figureHeight],'visible','off');
            
            h = 1 + this.showScarti;
            [ax, P] = tight_subplot(h,1,[this.padding 0],this.hMargins,this.wMargins);

            if this.showScarti
                set(ax(1),'Position',P{1}.*[1,0.78,1,1.3]);
                set(ax(2),'Position',P{2}.*[1,1,1,0.7]);
            end                   
                      
            set(ax(1),'box','on');
            hold(ax(1),'on');
            
            scatter(ax(1),this.datax, this.datay, "MarkerEdgeColor", this.dataColor);
                                       
            if isLinearFit
                if this.showModel
                    fplot(ax(1),@(x) this.par(1) + x * this.par(2), 'Color', this.modelColor, 'LineStyle', this.modelLineStyle, "LineWidth",this.lineWidth);
                end
                if this.showInitialParModel
                    
                    fplot(ax(1),@(x) this.previousPar(1) + x * this.previousPar(2), 'Color', 'k', 'LineStyle', this.modelLineStyle,"LineWidth",this.lineWidth);
                end             
            else
                if this.showModel
                    fplot(ax(1),@(x) this.model(this.par, x), 'Color', this.modelColor, 'LineStyle', this.modelLineStyle,"LineWidth",this.lineWidth);
                end
                if this.showInitialParModel                   
                    fplot(ax(1),@(x) this.model(this.previousPar, x), 'Color', 'k', 'LineStyle', this.modelLineStyle,"LineWidth",this.lineWidth);
                end
            end

            delta_x = abs(max(this.datax) - min(this.datax));
            delta_y = abs(max(this.datay) - min(this.datay));
            
            if this.xlim(1) >= this.xlim(2)
                % Gestione caso asse logaritmico per evitare il passaggio
                % dei limiti attraverso lo 0
                if this.logX
                    xlim(ax(1),sort([min(this.datax) * 0.2 max(this.datax) * 1.2]));
                else     
                    xlim(ax(1),sort([min(this.datax) - 0.1 * delta_x max(this.datax) + 0.1 * delta_x]));
                end
            else
                xlim(ax(1),[this.xlim(1) this.xlim(2)]);
            end

            if this.ylim(1) >= this.ylim(2)  
                % Gestione caso asse logaritmico per evitare il passaggio
                % dei limiti attraverso lo 0
                if this.logY
                    ylim(ax(1),sort([min(this.datay) * 0.2 max(this.datay) * 1.2]));
                else     
                    ylim(ax(1),sort([min(this.datay) - 0.1 * delta_y max(this.datay) + 0.1 * delta_y]));
                end             
            else
                ylim(ax(1),[this.ylim(1) this.ylim(2)]);
            end
            
            title(ax(1),this.name);
            ylabel(ax(1),this.labely);

            if this.logX
                set(ax(1), 'XScale', 'log')
            end

            if this.logY
                set(ax(1), 'YScale', 'log')
            end

            if this.showGrid
                grid(ax(1),'minor');
            end

            if ~this.showScarti
                xlabel(ax(1),this.labelx);
            else
                set(ax(1), 'XTickLabel', []);
            end
                        
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

            if this.showScarti                
                
                set(ax(2),'box','on');
                hold(ax(2),'on');

                if isLinearFit
                    sigmaScarti = sqrt(this.sigmay .^ 2 + (this.par(2) * this.sigmax) .^ 2);
                else
                    % Derivata numerica del modello per propagazione incertezze                   
                    sigma_propagata = propagation1D(@(x) this.model(this.par, x), this.datax, this.sigmax);
                    sigmaScarti = sqrt(this.sigmay .^2 + (sigma_propagata).^2);    
                end
                
                scarto_y = this.datay - this.yfit;
                               
                if this.logX
                    set(ax(2), 'XScale', 'log')
                end

                if this.showGrid
                    grid(ax(2),'minor');
                end

                % Linea orizzontale lungo y=0 negli scarti
                % Gestione caso asse logaritmico per evitare il passaggio
                % delle coordinate attraverso lasse 0
                if this.logX                   
                    x = sort([min(this.datax) * 0.2 max(this.datax) * 1.2]);
                    y = [0 0];
                else
                    x = [min(this.datax) - 0.1 * delta_x max(this.datax) + 0.1 * delta_x];
                    y = [0 0];
                end                
                line(ax(2),x, y, 'Color', this.modelColor, 'LineStyle', '-',"LineWidth", this.lineWidth)
                
                if this.continuosData
                    plot(ax(2),this.datax, this.datay,"LineWidth",this.lineWidth,"Color",[this.dataColor 1]);
                    plot(ax(2),this.datax, this.datay - sigmaScarti, "LineStyle","--","Color", [this.dataColor .3],"LineWidth",this.lineWidth);
                    plot(ax(2),this.datax, this.datay + sigmaScarti, "LineStyle","--","Color", [this.dataColor .3],"LineWidth",this.lineWidth);
                else
                    e2 = errorbar(ax(2),this.datax, scarto_y, sigmaScarti);
                    e2.LineStyle = 'none';                  
                    scatter(ax(2),this.datax, scarto_y, "MarkerEdgeColor", this.dataColor);
                end
               
                if this.xlim(1) >= this.xlim(2)
                    % Gestione caso asse logaritmico per evitare passaggio
                    % dei limiti attraverso lo 0
                    if this.logX
                        xlim(ax(2),sort([min(this.datax) * 0.2 max(this.datax) * 1.2]));
                    else     
                        xlim(ax(2),sort([min(this.datax) - 0.1 * delta_x max(this.datax) + 0.1 * delta_x]));
                    end
                else
                    xlim(ax(2),[this.xlim(1) this.xlim(2)]);
                end
                
                if this.resylim(1) >= this.resylim(2)
                    tmp_limits = max([abs(max(scarto_y)) + 2*abs(max(sigmaScarti)) abs(min(scarto_y)) + 2*abs(max(sigmaScarti))]);
                    ylim(ax(2),[-tmp_limits tmp_limits]); % Asse Y scarti simmetrico rispetto allo 0      
                else
                    ylim(ax(2),[this.resylim(1) this.resylim(2)])
                end
                
                ylabel(ax(2),this.reslabely);
                xlabel(ax(2),this.labelx);                                              
            end

            % Grafico errore zoom  ----------------------------------------

            if (this.showZoom)
                id = round(length(this.datax) / 2);
                x = this.datax(id);
                y = this.datay(id);
                ax(length(ax)+1) = axes("Position", this.zoomPosition);
                errorbar(ax(length(ax)),x, y, -this.sigmay(id), this.sigmay(id), -this.sigmax(id), this.sigmax(id), 'o');
                %xlim(ax(length(ax)), sort([(x - this.sigmax(id) * 1.5) (x + this.sigmax(id) * 1.5)]));
                %ylim(ax(length(ax)), sort([(y - this.sigmay(id) * 1.5) (y + this.sigmay(id) * 1.5)]));
                if this.showGrid
                    grid(ax(length(ax)),'minor');
                end
            end

            for a = ax
                set(a, "FontSize", this.fontSize);
            end

        end        
    end
end
