classdef fourierTransform < handle
    % ---------------------------------------------------
    % GITHUB: https://github.com/Sparisi01/matlab_lab_utility
    % ---------------------------------------------------
    % DEPENDENCIES:
    % - ./utils/exportFigure.m
    % - ./functionFit.m
    % ---------------------------------------------------

    properties
        
        data (:,1) double {mustBeReal, mustBeFinite} 
        sigmaData (:,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} 
        dt (1,1) double {mustBeReal, mustBeFinite}
        tollerance (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative}  % Tolleranza sulle ampiezze per evitare errori numerici. Porre a inf per disattivare.
        verbose (1, :) logical

        % Parametri che vengono riempiti dopo aver chiamato transform con i risultati        
        frequencies (:,1) double {mustBeReal, mustBeFinite}      
        amps (:,1) double {mustBeReal, mustBeFinite}
        phases (:,1) double {mustBeReal, mustBeFinite}
        dF (1,1) double {mustBeReal, mustBeFinite}
        sigmaAmps (:,1) double
        sigmaPhases (:,1) double
        
        % Peak detection
        peak_detection_centro_index (1,1) double 
        peak_detection_interval_index (1,1) double 
        peak_mean (1,1) double
        peak_sigma (1,1) double
        peak_amp (1,1) double
        peak_phase (1,1) double
        peak_amp_sigma (1,1) double
        peak_phase_sigma (1,1) double

        % Aesthetics
        name (1,1) string
        xLabel (1,1) string
        yLabel_phase (1,1) string
        yLabel_abs (1,1) string
        color (1,3) double {mustBeReal, mustBeFinite}
        xAxisLim (1,2) double {mustBeReal, mustBeFinite}
        xAxisAsOmegas (1,1) logical      
        mirrored (1,1) logical      
        fontSize (1, 1) double {mustBeReal, mustBeFinite}
        figureWidth (1, 1) double {mustBeReal, mustBeFinite}
        figureHeight (1, 1) double {mustBeReal, mustBeFinite}  
        stemPlot (1,1) logical  
    end

    methods
        
        % -----------------------------------------------------------------
        
        function this = fourierTransform()
            this.data = []; 
            this.sigmaData = 0;
            this.frequencies = [];
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
            this.stemPlot = 1;
        end

        
        % https://www.gaussianwaves.com/2015/11/interpreting-fft-results-obtaining-magnitude-and-phase-information/           
        function [frequencies, amps, phases, sigmaAmps, sigmaPhases] = transform(this)
            arguments
                this                
            end
                                
            fft_data = fftshift(fft(this.data))/length(this.data);                                  
            Fs = 1/this.dt;                                            
            this.dF = Fs/length(this.data); 
            this.frequencies = (-Fs/2:this.dF:Fs/2-this.dF)';
            this.amps = abs(fft_data);            
            
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
            
            % Considero solo le sigma associate a x positivo e applico il cambio variabile x -> x/2;
            tmp_sigmaData = zeros(size(this.sigmaData));
            for ii = 1:(length(this.sigmaData)/2)
                tmp_sigmaData(ii) =  linearSampling(1:length(this.sigmaData),this.sigmaData,ii/2 + (length(this.sigmaData)/2));
            end

            FFT_var_prime = 0.5 * fftshift(fft(tmp_sigmaData.^2))/length(tmp_sigmaData);
            dF_prime = this.dF*2;
            
            FFT_var = linearSampling((-Fs:dF_prime:Fs-dF_prime)', FFT_var_prime, this.frequencies);

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
                showFig (1,1) logical = 1
            end
            
            [frequencies, amps, phases, sigmaAmps, sigmaPhases] = this.transform();

            [fig, ax] = plotTransform(this, amps, sigmaAmps,1);
                        
            if showFig
                set(fig, 'visible', 'on'); 
            end

            if (strlength(fileName) > 0)
                exportFigure(fig, ax, fileName,this.fontSize, this.figureWidth, this.figureHeight);
            end 
        end

        % -----------------------------------------------------------------
        
        function [frequencies, phases, amps, sigmaAmps, sigmaPhases, fig, ax] = plotPhaseTransform(this, fileName, showFig)
            arguments
                this,
                fileName (1,1) string = "",
                showFig (1,1) logical = 1
            end
                        
            [frequencies, amps, phases, sigmaAmps, sigmaPhases] = this.transform();
            [fig, ax] = plotTransform(this, phases, sigmaPhases, 0);
            
            if showFig
                set(fig, 'visible', 'on'); 
            end

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
            
            [~,index_zero] = min(abs(this.frequencies));
                      
            ampiezza_meta_altezza = 0;

            if this.peak_detection_centro_index == inf 
                 [MAX,I] = max(this.amps(index_zero:end));
                 index_centro = index_zero + I;
                 ampiezza_meta_altezza = MAX/2;
            end

            if this.peak_detection_interval_index == inf    
                intervallo = 5;                      
            end
            
            intervallo_up = index_centro + intervallo;
            intervallo_do = index_centro - intervallo;     
                           
            
            tmp_freq = this.frequencies(intervallo_do:intervallo_up);
            tmp_amps = this.amps(intervallo_do:intervallo_up);
            % tmp_s_amps = this.sigmaAmps(intervallo_do:intervallo_up);
            tmp_phases = this.phases(intervallo_do:intervallo_up);

            x_meta_altezza_up = linearSampling(this.amps(index_centro-1:index_centro+intervallo),this.frequencies(index_centro-1:index_centro+intervallo),ampiezza_meta_altezza);
            x_meta_altezza_down = linearSampling(this.amps(index_centro-intervallo:index_centro-1),this.frequencies(index_centro-intervallo:index_centro-1),ampiezza_meta_altezza);
            
            this.peak_mean = (x_meta_altezza_up + x_meta_altezza_down)/2;
            this.peak_sigma = this.dF/sqrt(12) * sqrt(2);

            peak_mean = this.peak_mean;
            peak_sigma = this.peak_sigma;

            % Interpolazione lineare ampiezza picco
            peak_amp = linearSampling(tmp_freq,tmp_amps,peak_mean);
            peak_amp_sigma = this.sigmaAmps(index_centro);
            
            % Stima fase con nearest point
            peak_phase = tmp_phases(intervallo+1);
            peak_phase_sigma = this.sigmaPhases(index_centro);
            
            this.peak_amp = peak_amp;
            this.peak_amp_sigma = peak_amp_sigma;
            this.peak_phase = peak_phase;
            this.peak_phase_sigma = peak_phase_sigma;
            
        end        
    end

    methods (Hidden)
        function [fig, ax] = plotTransform(this, yData, sigmaY, isAbsPlot)
            arguments
                this,               
                yData (:,1) double {mustBeReal, mustBeFinite}, 
                sigmaY (:,1) double {mustBeReal, mustBeFinite}, 
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
            
            
            if this.stemPlot
                stem(xAxisData,yData,"filled", "Color", this.color);           
            else
                errorbar(xAxisData,yData,sigmaY,"Color",this.color,"Marker",".","LineStyle", "none");
            end
            
            
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
