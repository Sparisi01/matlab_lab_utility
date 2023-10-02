classdef functionFit < handle

    properties
        datax (:, 1) double {mustBeReal, mustBeFinite}
        datay (:, 1) double {mustBeReal, mustBeFinite}
        sigmax (:, 1) double {mustBeReal, mustBeFinite}
        sigmay (:, 1) double {mustBeReal, mustBeFinite}
        model (1, 1)
        par (:, 1) double {mustBeReal, mustBeFinite}
        errpar (:, 1) double
        ub (:, 1) double
        lb (:, 1) double
        units (:, 1) string
        parnames (:, 1) string
        yfit (:,1) double {mustBeReal, mustBeFinite}
        chi2norm (1,1) double {mustBeNonnegative}
        fig (1,1)  
        axes (:,1) 
        dof (1,1) double {mustBeNonnegative}
        pValue (1,1) double {mustBeNonnegative}
        toShowPar (:,1) logical
        showChi (1,1) logical
        showBox (1,1) logical 
        showScarti (1, 1) logical
        showGrid (1, 1) logical
        showZoom (1, 1) logical      
        zoomPosition (1, 4) double {mustBeReal, mustBeFinite} 
        modelColor (1, 3) double {mustBeReal, mustBeFinite}
        dataColor (1, 3) double {mustBeReal, mustBeFinite}
        name (1, 1) string
        labelx (1, 1) string
        labely (1, 1) string
        logX (1,1) logical
        logY (1,1) logical
        boxPosition (1, 2) double {mustBeReal, mustBeFinite}
        pedice (1, 1) char    
        fontSize (1, 1) double {mustBeReal, mustBeFinite}
        figureWidth (1, 1) double {mustBeReal, mustBeFinite}
        figureHeight (1, 1) double {mustBeReal, mustBeFinite}
        nolog (1, :) logical
        verbose (1, :) logical
    end

    methods
        
        function self = functionFit()
            % Dati --------------------------
            self.datax = [];
            self.datay = [];
            self.sigmax = [];
            self.sigmay = [];
            % Parametri modello -------------
            self.model = @(par, x) par(1) + x * par(2); % Modello su cui fittare
            self.par = []; % Valore dei parametri
            self.errpar = []; % Errore dei parametri
            self.ub = []; % UpperBound parametri
            self.lb = []; % LowerBound parametri
            self.units = []; % Units for parameters in legend box
            self.parnames = []; % Parameters name in legend box
            % Risultati fit -----------------
            self.yfit = [];
            self.chi2norm = inf;
            self.dof = 0; 
            self.pValue = 0; 
            self.fig = figure(); % Fig dopo aver generato figura e residui
            self.axes = []; % Array contenenti gli assi dopo aver generato figura e residui
            % Estetica ----------------------
            self.toShowPar = []; % Scegli quali parametri mostrate e quali no con un array di bool della stessa dimensione di par
            self.showChi = 1; % Mostra o no chi quadro in legenda
            self.showBox = 1; % Mostra o no box parametri
            self.showScarti = 1; % Mostra o no grafico degli scarti
            self.modelColor = [1 0 0];
            self.dataColor = [0.00 0.45 0.74];
            self.nolog = false;
            self.verbose = false;
            self.name = "Model"; % Titolo grafico
            self.labelx = "X Axes";
            self.labely = "Y Axes";        
            self.logX = 0;
            self.logY = 0;
            self.boxPosition = [0.55, 0.55]; % [x, y]
            self.pedice = ' '; % Pedice parametri legenda. Utile se si hanno molti grafici con parametri omonomi.
            self.showZoom = false; % Mostra grafico con zoom su un punto e barre incertezza
            self.zoomPosition = [0.21, 0.75, 0.15, 0.15]; % [x, y, w, h]
            self.showGrid = 1;
            self.fontSize = 14;
            self.figureWidth = 8; % Larghezza immagine salvata in pollici
            self.figureHeight = 6; % Altezza immagine salvata in pollici
        end
        
        % Genera immagine plot usando il modello dato
        function [par, errpar, yfit, chi2norm, dof, pValue, fig] = plotModelFit(self, file_name)

            arguments
                self
                file_name (1, 1) string = ""
            end

            % Esegui fit lineare e salva negli attributi dell'oggetto i nuovi valori dei parametri
            [par, errpar, yfit, chi2norm, dof, pValue] = modelFit(self);
            self.par = par;
            self.errpar = errpar;
            self.chi2norm = chi2norm;
            self.yfit = yfit;
            self.dof = dof;
            self.pValue = pValue;

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

            % Esegui fit non lineare e salva negli attributi dell'oggetto i nuovi valori dei parametri
            [par, errpar, yfit, chi2norm, dof, pValue] = linearFit(self);
            self.par = par;
            self.errpar = errpar;
            self.chi2norm = chi2norm;
            self.yfit = yfit;
            self.dof = dof;
            self.yfit = yfit;
            self.dof = dof;
            self.pValue = pValue;

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
            if ~safetyCheck(self)
                return
            end

            data_x = self.datax;
            data_y = self.datay;

            % Inizializza le sigma unitarie se non impostate
            if (isempty(self.sigmax))
                sigma_x = ones(size(self.datax));
                warning("Incertezze su datax non assegnate.")
            else
                sigma_x = self.sigmax;
            end   

            if (isempty(self.sigmay))            
                sigma_y = ones(size(self.datay));
                warning("Incertezze su datay non assegnate.")
            else
                sigma_y = self.sigmay;
            end      
           

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

            yfit = data_x * b + a;
            par(1) = a;
            par(2) = b;
            errpar(1) = sigma_a;
            errpar(2) = sigma_b;
            chi2norm = chi2 / (dataN - 2);
            dof = dataN - 2;
            pValue = chi2cdf(chi2, dof);
        end

        % Funzione per fit non lineare
        function [par, errpar, yfit, chi2norm, dof, pValue] = modelFit(self)

            % Contolla validità parametri
            if ~safetyCheck(self)
                return
            end

            options = optimset('lsqnonlin');

            % Definizione funzione scarti a partire dal modello
            scarti = @(par, xd, yd, ed) (self.model(par, xd) - yd) ./ ed;

            % Aggiusta dimensioni upper e lower bound. Se non impostate
            % vengono inizializzati a due vettori delle dimensione di par a
            % valori + e - inf
            if (isempty(self.ub))
                tmp_ub = ones(size(self.par)) * inf;
            else
                tmp_ub = self.ub;
            end

            if (isempty(self.lb))
                tmp_lb = ones(size(self.par)) * -inf;
            else
                tmp_lb = self.lb;
            end         
            
            % Inizializza le sigma unitarie se non impostate
            if (isempty(self.sigmax))
                tmp_sigmax = ones(size(self.datax));
                warning("Incertezze su datax non assegnate.")
            else
                tmp_sigmax = self.sigmax;
            end   

            if (isempty(self.sigmay))
                tmp_sigmay = ones(size(self.datay));
                warning("Incertezze su datay non assegnate.")
            else
                tmp_sigmay = self.sigmay;
            end      

            [par, resnorm, ~, ~, ~, ~, jacobian] = lsqnonlin(scarti, self.par, tmp_lb, tmp_ub, options, self.datax, self.datay, tmp_sigmay);

            %Covariance Matrix
            covar = inv(jacobian' * jacobian);
            %Variance
            var = diag(covar);
            sigma = sqrt(var);
            sigmaf = full(sigma);
            dof = (length(self.datax) - length(par));
            chi2norm = resnorm / dof;
            errpar = sigmaf * sqrt(chi2norm);
            yfit = self.model(par, self.datax);
            pValue = chi2cdf(resnorm, dof);
        end
       
    end

    methods (Hidden)

        % Fa tante cose belle
        function check = safetyCheck(self)           
            check = true;
            % Controlla dimensione dati
            if any(size(self.datay) ~= size(self.datax)) || ...
                    any(size(self.datay) ~= size(self.sigmay)) || ...
                    any(size(self.sigmax) ~= size(self.sigmay))
                check = false;
                error( 'u:stuffed:it' , 'datax, datay, sigmax. and sigmay devono essere della stesa dimensione'); 
            end

            % Controlla dimensione bound
            if ~(isempty(self.ub) || isempty(self.lb)) && ...
                    any(size(self.ub) ~= size(self.lb))
                check = false;
                error( 'u:stuffed:it' , "ub e lb devono essere della stesa dimensione");                
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

            delta_x = max(self.datax) - min(self.datax);
            
            % Disegna la funzione con parametri ottimizzati.
            if isLinearFit
                fplot(@(x) self.par(1) + x * self.par(2), 'Color', self.modelColor, 'LineStyle', '-')
            else
                fplot(@(x) self.model(self.par, x), 'Color', self.modelColor, 'LineStyle', '-')
            end

            xlim([min(self.datax) - 0.1 * delta_x max(self.datax) + 0.1 * delta_x]);
            ylim([min(self.datay) + 0.1 * min(self.datay) max(self.datay) + 0.1 * max(self.datay)]);

            if self.showGrid
                    grid minor;
            end
            hold on;
            % Plotta in punti
            scatter(self.datax, self.datay, "MarkerEdgeColor", self.dataColor);
            
            title(self.name);           
            ylabel(self.labely);

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
                    tmp_units = blanks(length(self.par));
                else
                    tmp_units = self.units;
                end
    
                % Aggiusta array nome parametri. Se non impostato viene inizializzato a un vettore delle
                % dimensioni di par formato da caratteri crescenti in ordine
                % alfabetico
                if (isempty(self.parnames))
                    tmp_parnames = char('a' + (1:length(self.par))-1);
                else                
                    tmp_parnames = self.parnames;
                end
                
                if isempty(self.toShowPar)               
                    self.toShowPar = ones(size(self.par));
                else 
                    if any(size(self.toShowPar) ~= size(self.par))
                        self.toShowPar = ones(size(self.par));
                        warning("Dimensione di par diversa da toShowPar.")
                    end
                end
                % Costruisci dinamicamente la legenda con n parametri.
                txt = "";          
                
                for ii = 1:length(self.par)
                    if self.toShowPar(ii)
                        t = numberToText(self.par(ii), self.errpar(ii));
        
                        if (self.pedice ~= ' ')
                            txt(length(txt) + 1) = tmp_parnames(ii) + "_{" + self.pedice + "} = " + t + " " + tmp_units(ii);
                        else
                            txt(length(txt) + 1) = tmp_parnames(ii) + " = " + t + " " + tmp_units(ii);
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
    
                e2 = errorbar(self.datax, scarto_y, sigmaScarti);
                e2.LineStyle = 'none';
                xlim([min(self.datax) - 0.1 * delta_x max(self.datax) + 0.1 * delta_x]);
                ylim([-max(scarto_y .* 2) - max(sigmaScarti) max(scarto_y .* 2) + max(sigmaScarti)]);
                x = [min(self.datax) - 0.1 * delta_x max(self.datax) + 0.1 * delta_x];
                y = [0 0];
                line(x, y, 'Color', self.modelColor, 'LineStyle', '-')
                if self.showGrid
                    grid minor;
                end
                hold on;
                scatter(self.datax, scarto_y, "MarkerEdgeColor", self.dataColor);
    
                % title("Residui da modello lineare", "fontSize", self.fontSize);
                ylabel("Residui");
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
                axes("Position", self.zoomPosition);
                errorbar(x, y, -self.sigmay(id), self.sigmay(id), -self.sigmax(id), self.sigmax(id), 'o');
                xlim([(x - self.sigmax(id) * 1.5) (x + self.sigmax(id) * 1.5)]);
                ylim([(y - self.sigmay(id) * 1.5) (y + self.sigmay(id) * 1.5)]);
                if self.showGrid
                    grid minor;
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
