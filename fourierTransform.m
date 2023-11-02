% -------------------------------------------------
% Funzione per eseguire agilmente la trasformata
% di fourier di un set di dati e generare il grafico
% in ampiezza e fase.
% ---------------------------------------------------
% DIPENDENZE:
% - ./utils/exportFigure.m
% - ./functionFit.m
% ---------------------------------------------------

% TODO
% Formalizzare propagazione incerteze

classdef fourierTransform < handle
    properties
        
        data (:,1) double {mustBeReal, mustBeFinite} 
        sigmaData (:,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} 
        times (:,1) double {mustBeReal, mustBeFinite} 
        dt (1,1) double {mustBeReal, mustBeFinite} % Intervallo di campionamento
        tollerance (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} % Tolleranza sulle ampiezze per evitare errori numerici. Porre a inf per disattivare.
        verbose (1, :) logical

        % Parametri che vengono riempiti dopo aver chiamato transform con i risultati        
        frequencies (:,1) double {mustBeReal, mustBeFinite} % Vettore frequenze      
        amps (:,1) double {mustBeReal, mustBeFinite} % Vettore ampiezze
        phases (:,1) double {mustBeReal, mustBeFinite} % Vettore fasi
        dF (1,1) double {mustBeReal, mustBeFinite} % intervallo di frequenza minimo risolvibile
        sigmaAmps (:,1) double % Vettore incertezze sulle ampiezze
        sigmaPhases (:,1) double % Vettore incertezze sulle fasi
        
        % Parametri peak detection
        peak_detection_centro_index (1,1) double 
        peak_detection_interval_index (1,1) double 
        peak_mean (1,1) double % Valore medio picco
        peak_sigma (1,1) double % Sigma picco
        peak_amp (1,1) double % Ampiezza picco
        peak_phase (1,1) double % Fase picco
        peak_amp_sigma (1,1) double
        peak_phase_sigma (1,1) double

        % Parametri estetici
        name (1,1) string
        xLabel (1,1) string
        yLabel_phase (1,1) string
        yLabel_abs (1,1) string
        color (1,3) double {mustBeReal, mustBeFinite}
        xAxisLim (1,2) double {mustBeReal, mustBeFinite}
        xAxisAsOmegas (1,1) logical % Visualizza pulsazioni sull'asse x invece che frequenze       
        mirrored (1,1) logical % Visualizza grafico specchiato       
        fontSize (1, 1) double {mustBeReal, mustBeFinite} % Dimensione font 
        figureWidth (1, 1) double {mustBeReal, mustBeFinite} % Larghezza immagine salvata in pollici
        figureHeight (1, 1) double {mustBeReal, mustBeFinite} % Altezza immagine salvata in pollici      
    end

    methods
        
        % -----------------------------------------------------------------
        
        function this = fourierTransform()
            this.data = []; 
            this.sigmaData = [];
            this.frequencies = [];
            this.times = [];
            this.dt = 1; 
            this.dF = 0;
            this.sigmaAmps = [];
            this.sigmaPhases = [];
            this.peak_detection_centro_index = inf;
            this.peak_detection_interval_index = inf;
            this.peak_mean = inf;
            this.peak_sigma = inf;
            this.peak_amp = inf;
            this.peak_amp_sigma = inf;
            this.peak_phase = inf;
            this.peak_phase_sigma = inf;
            this.name = "FFT";
            this.xLabel = "[Hz]";
            this.yLabel_phase = "\angle{FFT} [deg]";
            this.yLabel_abs = "|FFT| [V]";
            this.color = [0 0 1];
            this.mirrored = 0;
            this.xAxisLim = [0,0];
            this.xAxisAsOmegas = 0;
            this.tollerance = 1e-5;
            this.fontSize = 14; 
            this.figureWidth = 8; 
            this.figureHeight = 6; 
            this.verbose = 1;
        end

        % -----------------------------------------------------------------
        % Utilizzo della funzione FFT per estrarre vettore di ampiezze e fasi dal vettore dati
        % https://www.gaussianwaves.com/2015/11/interpreting-fft-results-obtaining-magnitude-and-phase-information/           
        function [frequencies, amps, phases, sigmaAmps, sigmaPhases] = transform(this)
            arguments
                this                
            end
                     
            % FFT dati
            fft_data = fftshift(fft(this.data))/length(this.data);
            % Frequenza massima risolvibile                      
            Fs = 1/this.dt;                                
            % Incertezza sulle frequenze
            this.dF = Fs/length(this.data); 
           
            % Vettore frequenze
            this.frequencies = (-Fs/2:this.dF:Fs/2-this.dF)';
                       
            % Vettore Ampiezze
            this.amps = abs(fft_data);            
            
            % Permetti impostare la sigma costante per tutto il set di dati
            % passando uno scalare
            if length(this.sigmaData) == 1
                this.sigmaData = ones(size(this.data)) * this.sigmaData;            
            end

            this.sigmaData = abs(this.sigmaData);

            % Funzione per gestire le divisioni 0/0
            % necessaria dopo aver inserito il treshold sulle ampiezze
            function z = custom_division(x,y)                 
                if (x == 0) + (y == 0) == 2 
                    z = 1;               
                else 
                    z = x./y; 
                end
            end

            % Calcolo incertezze tramite propagazione della serie di Fourier                                                
            % Vettore ausiliario incertezze sulle frequenze
            sigmaTimes = ones(size(this.times)) * this.dt/sqrt(12);

            %FFT_var_tmp = fftshift(fft(this.sigmaData.^2))/length(this.sigmaData);
            first_integral = this.dt * cumtrapz(this.sigmaData.^2.*sigmaTimes.^2);
            second_integral = this.dt * cumtrapz(this.times.^2.*this.sigmaData.^2);
            third_integral = this.dt * cumtrapz(this.data.^2.*sigmaTimes.^2);   
            forth_integral = this.dt * cumtrapz(this.sigmaData.^2); 
            correzione_primo_ordine = (this.frequencies/sqrt(2*pi))*(first_integral(end) + second_integral(end) + third_integral(end));
            correzione_primo_ordine = 0;
            FFT_var = forth_integral + correzione_primo_ordine;
            
            %FFT_var = fftshift(fft(this.sigmaData.^2))/length(this.sigmaData);
            FFT_sigma_real = sqrt(abs(real(FFT_var)));
            FFT_sigma_imag = sqrt(abs(imag(FFT_var)));
            FFT_real = real(fft_data);
            FFT_imag = imag(fft_data);       
            

            % Incertezza ampiezza e fase trasformata
            this.sigmaAmps = 1./(sqrt(FFT_real.^2 + FFT_imag.^2)) .* sqrt((FFT_real.*FFT_sigma_real).^2 + (FFT_imag.*FFT_sigma_imag).^2); 
            this.sigmaPhases = 1./(1+custom_division(FFT_imag,FFT_real).^2).*sqrt((1./FFT_real .* FFT_sigma_imag).^2 + (FFT_imag./FFT_real.^2 .* FFT_sigma_real).^2)*180/pi;
                   
            % Vettore Fasi con trasholds per rimuovere errori di computazione  
            tmp_fft_data = fft_data;   
            threshold = max(abs(fft_data))*this.tollerance;
            tmp_fft_data(abs(fft_data)<threshold) = 0;  
            this.phases = atan2(imag(tmp_fft_data),real(tmp_fft_data))*180/pi; 

            % Output
            amps = this.amps;
            phases = this.phases;
            frequencies = this.frequencies;  
            sigmaAmps = this.sigmaAmps;
            sigmaPhases = this.sigmaPhases;   
        end
        
        % -----------------------------------------------------------------
        
        function [frequencies, phases, amps,sigmaAmps,sigmaPhases, fig, ax] = plotAbsTransform(this, fileName, showFig)
            arguments
                this,
                fileName (1,1) string = "",
                showFig (1,1) logical = 0
            end
            
            [frequencies, amps, phases, sigmaAmps, sigmaPhases] = this.transform();

            [fig, ax] = plotTransform(this, amps,1);
                        
            % Rendi visibile figura
            if showFig
                set(fig, 'visible', 'on'); 
            end

            % Esporta figura
            if (strlength(fileName) > 0)
                exportFigure(fig, ax, fileName,this.fontSize, this.figureWidth, this.figureHeight);
            end 
        end

        % -----------------------------------------------------------------
        
        function [frequencies, phases, amps, sigmaAmps, sigmaPhases, fig, ax] = plotPhaseTransform(this, fileName, showFig)
            arguments
                this,
                fileName (1,1) string = "",
                showFig (1,1) logical = 0
            end
                        
            [frequencies, amps, phases, sigmaAmps, sigmaPhases] = this.transform();
            [fig, ax] = plotTransform(this, phases, 0);
            
            % Rendi visibile figura
            if showFig
                set(fig, 'visible', 'on'); 
            end

            % Esporta figura
            if (strlength(fileName) > 0)
                exportFigure(fig, ax, fileName,this.fontSize, this.figureWidth, this.figureHeight);
            end 
        end        

        % -----------------------------------------------------------------
        
        function [peak_mean, peak_sigma, peak_amp, peak_amp_sigma, peak_phase, peak_phase_sigma] = peakDetection(this)
            arguments
                this                              
            end
                        
            [~,~,~,~,~] = this.transform();
                   
            index_centro = this.peak_detection_centro_index;
            intervallo = this.peak_detection_interval_index;
            
            % Fissa valori di default se non impostati
            index_zero = round(length(this.amps)/2 + 1);
            if this.peak_detection_centro_index == inf 
                [~,I] = max(this.amps(index_zero:end));
                index_centro = index_zero + I;
            end

            if this.peak_detection_interval_index == inf    
                intervallo = 10;                      
            end
            
            % Definizione intervallo attorno al centro
            intervallo_up = index_centro + intervallo;
            intervallo_do = index_centro - intervallo;     
                           
            % Seleziona solo dati nell'intorno del picco selezionato 
            tmp_freq = this.frequencies(intervallo_do:intervallo_up);
            tmp_amps = this.amps(intervallo_do:intervallo_up);
            tmp_s_amps = this.sigmaAmps(intervallo_do:intervallo_up);
            tmp_phases = this.phases(intervallo_do:intervallo_up);
            
            % Calcolo media e dmedia per il picco       
            probability_density_freq = tmp_amps/sum(tmp_amps);
            s_probability_density_freq = (sum(tmp_amps) - tmp_amps)/(sum(tmp_amps)^2).*tmp_s_amps;
            peak_mean = sum(tmp_freq .* probability_density_freq);
            peak_sigma = sqrt(sum(s_probability_density_freq.^2) + this.dF^2/12); 
            
            this.peak_mean = peak_mean;
            this.peak_sigma = peak_sigma;

            % Interpolazione lineare ampiezza picco
            peak_amp = linearSampling(tmp_freq,tmp_amps,peak_mean);
            peak_amp_sigma = this.sigmaAmps(index_centro);
            
            % Stima fase con nearest point
            peak_phase = tmp_phases(intervallo);
            peak_phase_sigma = this.sigmaPhases(index_centro);
            
            this.peak_amp = peak_amp;
            this.peak_amp_sigma = peak_amp_sigma;
            this.peak_phase = peak_phase;
            this.peak_phase_sigma = peak_phase_sigma;
            
        end        
    end

    methods (Hidden)
        function [fig, ax] = plotTransform(this, yData, isAbsPlot)
            arguments
                this,               
                yData (:,1) double {mustBeReal, mustBeFinite},                
                isAbsPlot (1,1) logical = 1,
            end
                                               
            fig = figure('units','inch','position',[0,0,this.figureWidth,this.figureHeight],'visible','off');
            
            ax = axes();            
            hold on           
            axis padded
            box on
            grid minor           

            % Scegli se mostrare frequenze o pulsazioni
            if this.xAxisAsOmegas
                xAxisData = this.frequencies*2*pi;
                dX = this.dF*2*pi;                
            else
                xAxisData = this.frequencies;
                dX = this.dF;                                
            end
            xlabel(this.xLabel);
            
            % Plotta dati
            stem(xAxisData,yData,"filled", "Color", this.color);
            
            % Imposta limiti
            if this.xAxisLim == zeros(1,2) 
                if this.mirrored
                    xlim([min(xAxisData) - dX, max(xAxisData) + dX]);
                else
                    xlim([0, max(xAxisData) + dX]);
                end
            else
                xlim(sort(this.xAxisLim));
            end
            
            if isAbsPlot
                ylabel(this.yLabel_abs);
            else
                ylabel(this.yLabel_phase);
            end
            
            title(this.name);           
            ax.FontSize = 14;               
                                          
        end
    end
end
