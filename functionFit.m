classdef functionFit < handle
    % ---------------------------------------------------
    % GITHUB: https://github.com/Sparisi01/matlab_lab_utility
    % ---------------------------------------------------
    % DEPENDENCIES: 
    % - ./utils/numberToText.m
    % - ./utils/exportFigure.m
    % - ./utils/propagation1D.m
    % - ./utils/avoidOversampling.m
    % - ./utils/tight_subplot.m
    % - ./utils/textBox.m
    % - ./utils/geograficalPlacement.m
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
        yfit (:, 1) double {mustBeFinite}
        chi2norm (1, 1) double {mustBeNonnegative}
        fig (1, 1)
        axes (:, 1)
        dof double {mustBeNonnegative}
        pValue double {mustBeNonnegative}
        showParArray (:, 1) logical
        showChi logical
        showChiNorm logical
        showPValue logical
        showBox logical
        showRes logical
        showGrid logical
        showZoom logical
        showInitialParModel (1,1) logical 
        showModel logical
        continuosData logical
        zoomPosition string
        modelColor (1, 4) double {mustBeReal, mustBeFinite}
        modelLineStyle (1, 1) string
        dataColor (1, 3) double {mustBeReal, mustBeFinite}
        lineWidth double
        title string
        labelx string
        labely string
        reslabely (1, 1) string
        logX (1, 1) logical
        logY (1, 1) logical
        xlim (1, 2) double {mustBeReal, mustBeFinite}
        ylim (1, 2) double {mustBeReal, mustBeFinite}
        resylim (1, 2) double {mustBeReal, mustBeFinite}
        boxPosition string
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

        function this = functionFit()
            % Data --------------------------
            this.datax = [];
            this.datay = [];
            this.sigmax = [];
            this.sigmay = [];

            % Madel parameters -------------
            this.model = @(par, x) par(1) + x * par(2);
            this.par = [];
            this.errpar = [];
            this.upperBounds = [];
            this.lowerBounds = [];
            this.units = [];
            this.parnames = [];
            this.noOversampling = 0;

            % Fit results  ----------------
            % After a fit function is called these parameters are
            % filled with the current fit results.
            this.yfit = [];
            this.chi2norm = inf;
            this.dof = 0;
            this.pValue = 0;
            this.fig = 0;
            this.axes = [];
            this.previousPar = [];
            
            % Aesthetics -------------------
            this.showParArray = [1 1];
            this.showChi = 1;
            this.showChiNorm = 0;
            this.showPValue = 0;
            this.showBox = 1;
            this.showRes = 1;
            this.showInitialParModel = 0;
            this.showModel = 1;
            this.continuosData = 0;
            this.modelColor = [1 0 0 1];
            this.modelLineStyle = '-';
            this.dataColor = [0.00 0.00 1.00];
            this.lineWidth = 2;
            this.xlim = [0 0];
            this.ylim = [0 0];
            this.resylim = [0 0];
            this.verbose = 1;
            this.title = "Model";
            this.labelx = "X Axes";
            this.labely = "Y Axes";
            this.reslabely = "Scarti";
            this.logX = 0;
            this.logY = 0;
            this.boxPosition = "northwest";
            this.pedice = ' ';
            this.showZoom = false;
            this.zoomPosition = "northeast";
            this.showGrid = 1;
            this.fontSize = 14;
            this.figureWidth = 8;
            this.figureHeight = 6;    
            this.hMargins = [0.12 0.08];
            this.wMargins = [0.13 0.13]; 
            this.padding = 0.05; 
        end

        function [par, errpar, yfit, chi2norm, dof, pValue, fig, ax] = plotModelFit(this, file_name, showFig)
            arguments
                this
                file_name string = "",
                showFig logical = 1
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
                file_name string = "",
                showFig logical = 1
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
                [this.datax, this.datay, this.sigmax, this.sigmay] = avoidOversampling(this.datax, this.datay, this.sigmax, this.sigmay);
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
                % è minore dell'incertezza sul b
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
            
            safetyCheck(this)

            if(this.noOversampling)
                [this.datax, this.datay, this.sigmax, this.sigmay] = avoidOversampling(this.datax, this.datay, this.sigmax, this.sigmay);
            end           
            
            res_model = @(par, xd, yd, ed) (this.model(par, xd) - yd) ./ ed;

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
            
            [par, resnorm, ~, flag, ~, ~, jacobian] = lsqnonlin(res_model, this.par, tmp_lb, tmp_ub, options, this.datax, this.datay, this.sigmay);
                      
            dof = (length(this.datax) - length(par));  
            yfit = this.model(par, this.datax);
            
            s_x_prop = propagation1D(@(x) this.model(par, x), this.datax, this.sigmax);
            s_res = sqrt(this.sigmay .^2 + (s_x_prop).^2); 
            
            res_y = abs(this.datay - yfit);            
            chi2norm = sum((res_y./s_res).^2)/dof;                      
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

        function safetyCheck(this)            
            
            if any(size(this.datay) ~= size(this.datax))                   
                error('u:stuffed:it', "'datax' and 'datay' must be of the same dimension.");
            end
            
            if ~(isempty(this.upperBounds) || isempty(this.lowerBounds)) && ...
                    any(size(this.upperBounds) ~= size(this.lowerBounds))

                error('u:stuffed:it', "'ub' and 'lb' must be of the same dimension.");
            end
            
            if (isempty(this.sigmax)) % Check sigmax dimensions
                this.sigmax = zeros(size(this.datay));
                warning("'datax' not assigned. Zeros array used as default value.")           
            else
                if length(this.sigmax) == 1
                    this.sigmax = ones(size(this.datax)) * this.sigmax;
                end                
            end
            
            if (isempty(this.sigmay)) % Check sigmay dimensions
                this.sigmay = ones(size(this.datay));
                warning("'datax' not assigned. Unitary array used as default value.")     
            else
                if length(this.sigmay) == 1
                    this.sigmay = ones(size(this.datay)) * this.sigmay;
                end 
            end
        end

        function [fig,ax] = generatePlotFig(this, isLinearFit)

            fig = figure('units','inch','position',[0,0,this.figureWidth,this.figureHeight],'visible','off');
            
            h = 1 + this.showRes;
            [ax, P] = tight_subplot(h,1,[this.padding 0],this.hMargins,this.wMargins);

            if this.showRes
                set(ax(1),'Position',P{1}.*[1,0.78,1,1.3]);
                set(ax(2),'Position',P{2}.*[1,1,1,0.7]);
            end                   
                      
            set(ax(1),'box','on');
            hold(ax(1),'on');
            title(ax(1),this.title);
            ylabel(ax(1),this.labely);
                        
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
            
            % Manage x log axes to avoid passing througth zero
            if this.xlim(1) >= this.xlim(2)
                if this.logX
                    xlim(ax(1),sort([min(this.datax) * 0.2 max(this.datax) * 1.2]));
                else     
                    xlim(ax(1),sort([min(this.datax) - 0.1 * delta_x max(this.datax) + 0.1 * delta_x]));
                end
            else
                xlim(ax(1),[this.xlim(1) this.xlim(2)]);
            end

            % Manage y log axes to avoid passing througth zero
            if this.ylim(1) >= this.ylim(2)  
                if this.logY
                    ylim(ax(1),sort([min(this.datay) * 0.2 max(this.datay) * 1.2]));
                else     
                    ylim(ax(1),sort([min(this.datay) - 0.1 * delta_y max(this.datay) + 0.1 * delta_y]));
                end             
            else
                ylim(ax(1),[this.ylim(1) this.ylim(2)]);
            end
            

            if this.logX
                set(ax(1), 'XScale', 'log')
            end

            if this.logY
                set(ax(1), 'YScale', 'log')
            end

            if this.showGrid
                grid(ax(1),'minor');
            end

            if ~this.showRes
                xlabel(ax(1),this.labelx);
            else
                set(ax(1), 'XTickLabel', []);
            end
                        
            if this.showBox                             
                if (isempty(this.units))       
                    for ii = 1:length(this.par)
                        this.units(ii) = "";
                    end                  
                else
                    if any(size(this.units) ~= size(this.par))
                        for ii = 1:length(this.par)
                            this.units(ii) = "";
                        end
                        warning("The argument 'par' must have the same size of 'units'. Values corrected automatically.")
                    end
                end
                
                if (isempty(this.parnames))
                    for ii = 1:length(this.par)
                        this.parnames(ii) = 'a' + ii - 1;
                    end                
                else
                    if any(size(this.parnames) ~= size(this.par))
                        for ii = 1:length(this.par)
                            this.parnames(ii) = 'a' + ii - 1;
                        end
                        warning("The argument 'par' must have the same size of 'parnames'. Values corrected automatically.")
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

                % Build the textBox dynamically based on the number of
                % parameters and userpreferences on chi2
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
                
                % Round Chi2 
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
                                
                textBox(txt(2:end),this.boxPosition,ax(1),this.fontSize); % Place textBox
            end

            if this.showRes                
                
                set(ax(2),'box','on');
                hold(ax(2),'on');

                if isLinearFit
                    s_res = sqrt(this.sigmay .^ 2 + (this.par(2) * this.sigmax) .^ 2);
                else
                    % Numeric propagation of 'sigmax' through the optimized model                   
                    s_x_prop = propagation1D(@(x) this.model(this.par, x), this.datax, this.sigmax);
                    s_res = sqrt(this.sigmay .^2 + (s_x_prop).^2);    
                end
                
                res_y = this.datay - this.yfit;
                               
                if this.logX
                    set(ax(2), 'XScale', 'log')
                end

                if this.showGrid
                    grid(ax(2),'minor');
                end

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
                    plot(ax(2),this.datax, this.datay - s_res, "LineStyle","--","Color", [this.dataColor .3],"LineWidth",this.lineWidth);
                    plot(ax(2),this.datax, this.datay + s_res, "LineStyle","--","Color", [this.dataColor .3],"LineWidth",this.lineWidth);
                else
                    e2 = errorbar(ax(2),this.datax, res_y, s_res, "Color",this.dataColor);
                    e2.LineStyle = 'none';                  
                    scatter(ax(2),this.datax, res_y, "MarkerEdgeColor", this.dataColor);
                end
               
                % Manage x log axes to avoid passing througth zero
                if this.xlim(1) >= this.xlim(2)
                    if this.logX
                        xlim(ax(2),sort([min(this.datax) * 0.2 max(this.datax) * 1.2]));
                    else     
                        xlim(ax(2),sort([min(this.datax) - 0.1 * delta_x max(this.datax) + 0.1 * delta_x]));
                    end
                else
                    xlim(ax(2),[this.xlim(1) this.xlim(2)]);
                end
                
                if this.resylim(1) >= this.resylim(2)
                    tmp_limits = max([abs(max(res_y)) + 2*abs(max(s_res)) abs(min(res_y)) + 2*abs(max(s_res))]);
                    ylim(ax(2),[-tmp_limits tmp_limits]);
                else
                    ylim(ax(2),[this.resylim(1) this.resylim(2)])
                end
                
                ylabel(ax(2),this.reslabely);
                xlabel(ax(2),this.labelx);                                              
            end

            if (this.showZoom)
                id = round(length(this.datax) / 2);
                x = this.datax(id);
                y = this.datay(id);
                ax(length(ax)+1) = geograficPlacement(axes("Position", [0,0,0.15,0.15]), this.zoomPosition, ax(1));
                errorbar(ax(length(ax)),x, y, -this.sigmay(id), this.sigmay(id), -this.sigmax(id), this.sigmax(id), 'o');
                xlim(ax(length(ax)), sort([(x - this.sigmax(id) * 1.5) (x + this.sigmax(id) * 1.5)]));
                ylim(ax(length(ax)), sort([(y - this.sigmay(id) * 1.5) (y + this.sigmay(id) * 1.5)]));
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
