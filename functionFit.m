classdef functionFit < handle

    properties
        datax (:, 1) double {mustBeReal, mustBeFinite}
        datay (:, 1) double {mustBeReal, mustBeFinite}
        sigmax (:, 1) double {mustBeReal, mustBeFinite}
        sigmay (:, 1) double {mustBeReal, mustBeFinite}
        model (1, 1)
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
        zoomPosition (1, 4) double {mustBeReal, mustBeFinite}
        modelColor (1, 3) double {mustBeReal, mustBeFinite}
        modelLineStyle (1, 1) string
        dataColor (1, 3) double {mustBeReal, mustBeFinite}
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
    end

    methods

        % Valori di default per le opzioni
        function self = functionFit()
            % Dati --------------------------
            self.datax = [];
            self.datay = [];
            % Se le incertezze non vengono definite vengono inizializzate a
            % un vettore unitario, Se viene passato uno scalare la stessa
            % incertezza viene applicata a ogni punto.
            self.sigmax = [];
            self.sigmay = [];

            % Parametri modello -------------
            self.model = @(par, x) par(1) + x * par(2); % Modello su cui eseguire il fit
            self.par = []; % Valore dei parametri
            self.previousPar = []; % Valore dei parametri
            self.errpar = []; % Errore dei parametri
            self.upperBounds = []; % UpperBound parametri
            self.lowerBounds = []; % LowerBound parametri
            self.units = []; % Units for parameters in legend box
            self.parnames = []; % Parameters name in legend box
            
            % Risultati fit  ----------------
            self.yfit = [];
            self.chi2norm = inf;
            self.dof = 0;
            self.pValue = 0;
            self.fig = figure(); % Fig dopo aver generato figura e residui
            self.axes = []; % Array contenenti gli assi dopo aver generato figura e residui
            
            % Estetica ----------------------
            self.showParArray = [1 1]; % Scegli quali parametri mostrate con un array di bool della stessa dimensione di par
            self.showDataArray = []; % Scegli quali dati mostrate con un array di bool della stessa dimensione dei dati. I dati non mostrati in grafico contribuiscono comunque al fit
            self.showChi = 1; % Mostra o no chi quadro in legenda
            self.showChiNorm = 0; % Mostra chi2normalizzato
            self.showPValue = 0; % Mostra PValue
            self.showBox = 1; % Mostra o no box parametri
            self.showScarti = 1; % Mostra o no grafico degli scarti
            self.showInitialParModel = 0; % Mostra modello con parametri iniziali su grafico
            self.showModel = 1; % Mostra modello sul grafico
            self.modelColor = [1 0 0]; % Colore linea modello
            self.modelLineStyle = '-';
            %self.dataColor = [0.00 0.45 0.74]; % Colore dati
            self.dataColor = [0.00 0.00 1.00]; % Colore dati
            self.xlim = [0 0]; % Xlim, se uguali o in ordine sbagliato viene impostato in automatico
            self.ylim = [0 0]; % Ylim, se uguali o in ordine sbagliato viene impostato in automatico
            self.resylim = [0 0]; % Ylim residui, se uguali o in ordine sbagliato viene impostato in automatico
            self.verbose = false;
            self.name = "Model"; % Titolo grafico
            self.labelx = "X Axes";
            self.labely = "Y Axes";
            self.reslabely = "Scarti"; % Label y scarti
            self.logX = 0; % Asse X logaritmico (per entrambi i grafici)
            self.logY = 0; % Asse Y logaritmico
            self.boxPosition = [0.55, 0.55]; % [x, y] la dimensione di aggiusta in automatico
            self.pedice = ' '; % Pedice parametri legenda. Utile se si hanno molti grafici con parametri omonomi.
            self.showZoom = false; % Mostra grafico con zoom su un punto e barre incertezza
            self.zoomPosition = [0.21, 0.75, 0.15, 0.15]; % [x, y, w, h]
            self.showGrid = 1; % Mostra griglia minor sui grafici
            self.fontSize = 14; % Dimensione font sia nelle label che nella box
            self.figureWidth = 8; % Larghezza immagine salvata in pollici
            self.figureHeight = 6; % Altezza immagine salvata in pollici
        end

        % Genera immagine plot usando il modello dato
        function [par, errpar, yfit, chi2norm, dof, pValue, fig] = plotModelFit(self, file_name)

            arguments
                self
                file_name (1, 1) string = ""
            end
            
            self.previousPar = self.par;
            % Esegui fit lineare e salva negli attributi dell'oggetto i nuovi valori dei parametri
            [par, errpar, yfit, chi2norm, dof, pValue] = modelFit(self);
            
            

            % Genera Figura
            [fig, ax] = generatePlotFig(self, false);
            self.fig = fig;
            self.axes = ax;

            % Esporta figura in formato png
            if (strlength(file_name) > 0)
                exportFigure(fig, ax, file_name);
            end

        end

        % Genera immagine plot usando il modello lineare
        function [par, errpar, yfit, chi2norm, dof, pValue, fig] = plotLinearFit(self, file_name)

            arguments
                self
                file_name (1, 1) string = ""
            end
            
            self.previousPar = self.par;
            % Esegui fit non lineare e salva negli attributi dell'oggetto i nuovi valori dei parametri
            [par, errpar, yfit, chi2norm, dof, pValue] = linearFit(self);
            
            
            
            % Genera Figura
            [fig, ax] = generatePlotFig(self, true);
            self.fig = fig;
            self.axes = ax;

            % Esporta figura in formato png
            if (strlength(file_name) > 0)
                exportFigure(fig, ax, file_name);
            end

        end

        % Funzione per fit lineare
        function [par, errpar, yfit, chi2norm, dof, pValue] = linearFit(self)

            % Contolla validità parametri
            safetyCheck(self)

            data_x = self.datax;
            data_y = self.datay;
            sigma_x = self.sigmax;
            sigma_y = self.sigmay;

            % Inizializza le sigma unitarie se non impostate
            %if (isempty(self.sigmax))
            %    sigma_x = max([0.01*max(data_x)*ones(size(data_x)), 0.01*data_x],[],2);
            %    warning("Incertezze su datax non assegnate. Inizializzate all'1%")
            %else
            %    sigma_x = self.sigmax;
            %end

            %if (isempty(self.sigmay))
            %    sigma_y = max([0.01*max(data_y)*ones(size(data_y)), 0.01*data_y],[],2);
            %    warning("Incertezze su datay non assegnate. Inizializzate all'1%")
            %else
            %    sigma_y = self.sigmay;
            %end

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
            
            self.par = par;
            self.errpar = errpar;
            self.chi2norm = chi2norm;
            self.yfit = yfit;
            self.dof = dof;
            self.pValue = pValue;

        end

        % Funzione per fit non lineare
        function [par, errpar, yfit, chi2norm, dof, pValue, flag] = modelFit(self)

            % Contolla validità parametri
            safetyCheck(self)

            % Definizione funzione scarti a partire dal modello
            scarti = @(par, xd, yd, ed) (self.model(par, xd) - yd) ./ ed;

            % Aggiusta dimensioni upper e lower bound. Se non impostate
            % vengono inizializzati a due vettori delle dimensione di par a
            % valori + e - inf
            if (isempty(self.upperBounds))
                tmp_ub = ones(size(self.par)) * inf;
            else
                tmp_ub = self.upperBounds;
            end

            if (isempty(self.lowerBounds))
                tmp_lb = ones(size(self.par)) * -inf;
            else
                tmp_lb = self.lowerBounds;
            end
            
            %if (isempty(self.sigmax))
            %    tmp_sigmax = max([0.01*max(self.datax)*ones(size(self.datax)), 0.01*self.datax],[],2);
            %    warning("Incertezze su datax non assegnate. Inizializzate all'1%")
            %else
            %    if size(self.sigmay) == [1,1]
            %        self.sigmax = ones(size(self.datax)) * self.sigmax;
            %    end
            %    tmp_sigmax = self.sigmax;
            %end

            % if (isempty(self.sigmay))
            %     tmp_sigmay = max([0.01*max(self.datay)*ones(size(self.datay)), 0.01*self.datay],[],2);
            %     warning("Incertezze su datay non assegnate. Inizializzate all'1%")
            % else
            %     if size(self.sigmay) == [1,1]
            %         self.sigmay = ones(size(self.datay)) * self.sigmay;
            %     end
            %     tmp_sigmay = self.sigmay;               
            % end

            if self.verbose
                options = optimoptions('lsqnonlin');
            else
                 options = optimoptions('lsqnonlin', 'Display', 'none');
            end

            [par, resnorm, ~, flag, ~, ~, jacobian] = lsqnonlin(scarti, self.par, tmp_lb, tmp_ub, options, self.datax, self.datay, self.sigmay);
            
            
            %Covariance Matrix
            covar = inv(jacobian' * jacobian);
            %Variance
            var = diag(covar);
            sigma = sqrt(var);
            sigmaf = full(sigma);

            % Results
            dof = (length(self.datax) - length(par));
            chi2norm = resnorm / dof;
            errpar = sigmaf * sqrt(chi2norm);
            yfit = self.model(par, self.datax);
            pValue = chi2cdf(resnorm, dof);    

            self.par = par;
            self.errpar = errpar;
            self.chi2norm = chi2norm;
            self.yfit = yfit;
            self.dof = dof;
            self.pValue = pValue;
        end

    end

    methods (Hidden)

        % Fa tante cose belle
        function safetyCheck(self)
            
            % Controlla dimensione dati
            if any(size(self.datay) ~= size(self.datax))                   
                error('u:stuffed:it', 'datax, datay devono essere della stesa dimensione');
            end

            % Controlla dimensione bound
            if ~(isempty(self.upperBounds) || isempty(self.lowerBounds)) && ...
                    any(size(self.upperBounds) ~= size(self.lowerBounds))

                error('u:stuffed:it', "ub e lb devono essere della stesa dimensione");
            end
            
            % Controlla dimensione incertezze
            if (isempty(self.sigmax))
                self.sigmax = max([0.01*max(self.datax)*ones(size(self.datax)), 0.01*self.datax],[],2);
                warning("Incertezze su datax non assegnate. Inizializzate all'1%")           
            else
                if size(self.sigmax) == [1,1]
                    self.sigmax = ones(size(self.datax)) * self.sigmax;
                end                
            end

            if (isempty(self.sigmay))
                self.sigmay = max([0.01*max(self.datay)*ones(size(self.datay)), 0.01*self.datay],[],2);
                warning("Incertezze su datay non assegnate. Inizializzate all'1%")        
            else
                if size(self.sigmay) == [1,1]
                    self.sigmay = ones(size(self.datay)) * self.sigmay;
                end 
            end
        end

        function [fig, ax] = generatePlotFig(self, isLinearFit)

            fig = figure();

            if self.showScarti
                h = 3;
                tiledlayout(h, 1);

                % Costruzione grafico principale
                ax(1) = nexttile(1, [h - 1, 1]);
            else
                ax(1) = axes();
            end

            delta_x = abs(max(self.datax) - min(self.datax));
            hold on
            box on

            if (isempty(self.showDataArray))
                self.showDataArray = ones(size(self.datax));
            else
                if size(self.showDataArray) ~= size(self.datax)
                    error('u:stuffed:it', "showDataArray e datax, datay devono essere della stessa dimensione.");
                end
            end

            scatter(self.datax.*self.showDataArray, self.datay.*self.showDataArray, "MarkerEdgeColor", self.dataColor);
            
            
            % Disegna la funzione con parametri ottimizzati e parametri iniziali     
            if isLinearFit
                if self.showModel
                    fplot(@(x) self.par(1) + x * self.par(2), 'Color', self.modelColor, 'LineStyle', self.modelLineStyle);
                end
                if self.showInitialParModel
                    hold on;
                    fplot(@(x) self.previousPar(1) + x * self.previousPar(2), 'Color', 'k', 'LineStyle', self.modelLineStyle);
                end             
            else
                if self.showModel
                    fplot(@(x) self.model(self.par, x), 'Color', self.modelColor, 'LineStyle', self.modelLineStyle);
                end
                if self.showInitialParModel
                    hold on;
                    fplot(@(x) self.model(self.previousPar, x), 'Color', 'k', 'LineStyle', self.modelLineStyle);
                end
            end

            % Imposta dinamicamente i limiti
            if self.xlim(1) >= self.xlim(2)
                xlim([min(self.datax) - 0.1 * delta_x max(self.datax) + 0.1 * delta_x]);
            else
                xlim([self.xlim(1) self.xlim(2)]);
            end

            if self.ylim(1) >= self.ylim(2)
                ylim([min(self.datay) - 0.1 * abs(min(self.datay)) max(self.datay) + 0.1 *abs(max(self.datay))]);
            else
                ylim([self.ylim(1) self.ylim(2)]);
            end
            
            title(self.name);
            ylabel(self.labely);
            
            if self.showGrid
                grid minor;
            end 

            if self.logX
                set(ax(1), 'XScale', 'log')
            end

            if self.logY
                set(ax(1), 'YScale', 'log')
            end

            if ~self.showScarti
                xlabel(self.labelx);
            else
                set(gca, 'XTickLabel', []);
            end
                        
            if self.showBox
                
                
                % Aggiusta array unitÃ  di misura. Se non impostate le unitÃ  di misura
                % sono inizializzate a un vettore delle dimensioni di par
                % formato da numeri crescenti da 1
                if (isempty(self.units))       
                    for ii = 1:length(self.par)
                        self.units(ii) = "";
                    end
                    
                else
                    if any(size(self.units) ~= size(self.par))
                        for ii = 1:length(self.par)
                            self.units(ii) = "";
                        end
                        warning("Dimensione di par diversa da units. Parametro inizializzato in automatico.")
                    end
                end
                

                % Aggiusta array nome parametri. Se non impostato viene inizializzato a un vettore delle
                % dimensioni di par formato da caratteri crescenti in ordine
                % alfabetico
                if (isempty(self.parnames))
                    for ii = 1:length(self.par)
                        self.parnames(ii) = 'a' + ii - 1;
                    end                
                else
                    if any(size(self.parnames) ~= size(self.par))
                        for ii = 1:length(self.par)
                            self.parnames(ii) = 'a' + ii - 1;
                        end
                        warning("Dimensione di par diversa da parnames. Parametro inizializzato in automatico.")
                    end
                end
                

                if isempty(self.showParArray)
                    self.showParArray = ones(size(self.par));
                else
                    if any(size(self.showParArray) ~= size(self.par))
                        self.showParArray = ones(size(self.par));
                        warning("Dimensione di par diversa da showParArray. Parametro inizializzato in automatico.")
                    end
                end

                % Costruisci dinamicamente la legenda con n parametri.
                txt = "";            
                for ii = 1:length(self.par)
                    if self.showParArray(ii)
                        t = numberToText(self.par(ii), self.errpar(ii));
                        
                        if (self.pedice ~= ' ')
                            txt(length(txt) + 1) = self.parnames(ii) + "_{" + self.pedice + "} = " + t + " " + self.units(ii);
                        else
                            txt(length(txt) + 1) = self.parnames(ii) + " = " + t + " " + self.units(ii);
                        end
                    end
                end

                % Se il chi2 è minore di 2 vengono tenute 2 cifre significative
                % sennò si arrotonda all'intero.
                if self.chi2norm * self.dof > 2
                    nRound = 0;
                else
                    nRound = 2;
                end

                if self.showChi
                    txt(length(txt) + 1) = "\chi^2_{" + self.pedice + "} = " + round(self.chi2norm * self.dof, nRound) + "/" + self.dof;
                end

                if self.showChiNorm
                    txt(length(txt) + 1) = "\chi^2_{" + self.pedice + "} = " + round(self.chi2norm, 2);
                end

                if self.showPValue
                    txt(length(txt) + 1) = "pValue_{" + self.pedice + "} = " + round(self.pValue, 2);
                end

                % Dinamic position
                annotation("textbox", [self.boxPosition 0 0], ...
                    "BackgroundColor", [1, 1, 1], ...
                    "fontSize", self.fontSize, ...
                    "String", txt(2:end), ...
                    'FitBoxToText', 'on' ...
                );
            end

            if self.showScarti
                % Costruzione grafico degli scarti
                ax(2) = nexttile([1 1]);

                scarto_y = self.datay - self.yfit;

                % Costruisci incertezza totale proiettando le incertezze su x
                % attraverso le derivate del modello.
                if isLinearFit
                    sigmaScarti = sqrt(self.sigmay .^ 2 + (self.par(2) * self.sigmax) .^ 2);
                else
                    sigmaScarti = self.sigmay;
                end

                e2 = errorbar(self.datax.*self.showDataArray, scarto_y.*self.showDataArray, sigmaScarti);
                e2.LineStyle = 'none';

                % Imposta dinamicamente i limiti
                if self.xlim(1) >= self.xlim(2)
                    xlim([min(self.datax) - 0.1 * delta_x max(self.datax) + 0.1 * delta_x]);
                else
                    xlim([self.xlim(1) self.xlim(2)]);
                end

                if self.resylim(1) >= self.resylim(2)
                    ylim([-max(scarto_y .* 2) - max(sigmaScarti) max(scarto_y .* 2) + max(sigmaScarti)]);
                else
                    ylim([self.resylim(1) self.resylim(2)])
                end

                x = [min(self.datax) - 0.1 * delta_x max(self.datax) + 0.1 * delta_x];
                y = [0 0];
                line(x, y, 'Color', self.modelColor, 'LineStyle', '-')

                if self.showGrid
                    grid minor;
                end

                hold on;
                scatter(self.datax.*self.showDataArray, scarto_y.*self.showDataArray, "MarkerEdgeColor", self.dataColor);

                ylabel(self.reslabely);
                xlabel(self.labelx);

                if self.logX
                    set(ax(2), 'XScale', 'log')
                end

            end

            % Costruisci grafico zoom errori su un punto
            % Sono mostrati sia gli errori su x che y
            if (self.showZoom)
                hold on;
                id = round(length(self.datax) / 2);
                x = self.datax(id);
                y = self.datay(id);
                ax(length(ax)+1) = axes("Position", self.zoomPosition);
                errorbar(x, y, -self.sigmay(id), self.sigmay(id), -self.sigmax(id), self.sigmax(id), 'o');
                xlim(ax(length(ax)), [(x - self.sigmax(id) * 1.5) (x + self.sigmax(id) * 1.5)]);
                ylim(ax(length(ax)), [(y - self.sigmay(id) * 1.5) (y + self.sigmay(id) * 1.5)]);

                if self.showGrid
                    grid on;
                end

            end

            % Imposta dimensione font per ogni asse.
            % Richiesto perchè il tiled layout lavora con due assi
            % differenti
            for a = 1:length(ax)
                set(ax(a), "fontSize", self.fontSize);
            end

        end

        % Funzione export figura in formato png
        function exportFigure(self, image_figure, image_axes, image_name)

            % Imposta dimensione font per ogni asse
            for a = image_axes
                set(a, "fontSize", self.fontSize);
            end

            set(image_figure, 'PaperUnits', 'inches');
            set(image_figure, 'PaperSize', [self.figureWidth self.figureHight]);

            set(image_figure, 'InvertHardcopy', 'on');
            set(image_figure, 'PaperUnits', 'inches');
            set(image_figure, 'PaperPosition', [0, 0, self.figureWidth, self.figureHight]);

            % Salva file come PNG
            print(image_name, '-dpng', '-r300');

        end

        % Number to text
        function text = numberToText(x, sx)
            % x Ã¨ un numero
            % sx Ã¨ la sua incertezza
            % cifre Ã¨ il numero di cifre significative

            og = floor(log10(abs(x))); % ordine di grandezza
            % mantissa_x = round(x / (10^og) * 10^cifre) / (10^cifre);

            ogs = floor(log10(abs(sx)));
            sx = round(sx / (10 ^ ogs)) * (10 ^ ogs);
            cifre = strlength(string(sx / (10 ^ og))) - 2;

            mantissa_x = sprintf("%0." + cifre + "f", x / (10 ^ og));
            mantissa_sx = sx / (10 ^ og) + "";

            if og == 0
                text = "(" + mantissa_x + " \pm " + mantissa_sx + ")";
            else
                text = "(" + mantissa_x + " \pm " + mantissa_sx + ")\times{}10^{" + og + "}";
            end

        end

    end

end
