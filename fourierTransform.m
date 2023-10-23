% -------------------------------------------------
% Funzione per eseguire agilmente la trasformata
% di fourier di un set di dati e generare il grafico
% in ampiezza e fase.
% ---------------------------------------------------
% DIPENDENZE:
% - ./utils/exportFigure.m
% ---------------------------------------------------

% TODO
% Formalizzare propagazione incerteze

classdef fourierTransform < handle
    properties
        
        data (:,1) double {mustBeReal, mustBeFinite} 
        sigmaData (:,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} 
        dt (1,1) double {mustBeReal, mustBeFinite} % Intervallo di tempo tra i campionamenti x
        tollerance (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} % Tolleranza sulle ampiezze per evitare errori numerici. Porre a inf per disattivare.
        verbose (1, :) logical

        % Parametri che vengono riempiti dopo aver chiamato transform con i risultati        
        frequencies (:,1) double {mustBeReal, mustBeFinite} % Vettore frequenze      
        amps (:,1) double {mustBeReal, mustBeFinite} % Vettore ampiezze
        phases (:,1) double {mustBeReal, mustBeFinite} % Vettore fasi
        dF (1,1) double {mustBeReal, mustBeFinite} % intervallo di frequenza minimo risolvibile
        sigmaAmps (:,1) double % Vettore incertezze sulle ampiezze
        sigmaPhases (:,1) double % Vettore incertezze sulle fasi

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
            this.dt = 1; % Intervallo temporale tra i dati
            this.dF = 0;
            this.sigmaAmps = [];
            this.sigmaPhases = [];
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
            
            % Vettore Fasi   
            tmp_fft_data = fft_data;   
            threshold = max(abs(fft_data))*this.tollerance;
            tmp_fft_data(abs(fft_data)<threshold) = 0;  
            this.phases = atan2(imag(tmp_fft_data),real(tmp_fft_data))*180/pi; 
            
            % Calcolo incertezze tramite propagazione della serie di Fourier                                                
            
            % Permetti impostare la sigma costante per tutto il set di dati
            % passando uno scalare
            if length(this.sigmaData) == 1
                this.sigmaData = ones(size(this.data)) * this.sigmaData;            
            end

            % Funzione per gestire le divisioni 0/0
            % necessaria dopo aver inserito il treshold sulle ampiezze
            function z = custom_division(x,y)                 
                if (x == 0) + (y == 0) == 2 
                    z = 1;               
                else 
                    z = x./y; 
                end
            end
            
            % Incertezza ampiezze
            this.sigmaAmps = sqrt(2/(length(this.data)+1)) * mean(this.sigmaData) * ones(size(this.data));
            % Incertezza fasi
            this.sigmaPhases = (1/sqrt((2*length(this.data)+1)) * mean(this.sigmaData)) ./ ...
            (real(fft_data).*sqrt(1+custom_division(imag(fft_data),real(fft_data)).^2));
            
            % Output
            amps = this.amps;
            phases = this.phases;
            frequencies = this.frequencies;  
            sigmaAmps = this.sigmaAmps;
            sigmaPhases = this.sigmaPhases;   
        end
        
        % -----------------------------------------------------------------
        
        function [frequencies, phases, amps,sigmaAmps,sigmaPhases, fig, ax] = plotAbsTransform(this, fileName)
            arguments
                this,
                fileName (1,1) string = ""
            end
                        
            [frequencies, amps, phases, sigmaAmps, sigmaPhases] = this.transform();

            [fig, ax] = plotTransform(this, amps, fileName, 1);
        end

        % -----------------------------------------------------------------
        
        function [frequencies, phases, amps, sigmaAmps, sigmaPhases, fig, ax] = plotPhaseTransform(this, fileName)
            arguments
                this,
                fileName (1,1) string = ""
            end
                        
            [frequencies, amps, phases, sigmaAmps, sigmaPhases] = this.transform();
            [fig, ax] = plotTransform(this, phases, fileName, 0);
        end        

        % -----------------------------------------------------------------
        
        function [frequenze_media, frequenze_sigma] = peakIdentifier(this, intervallo)
            arguments
                this,
                intervallo (:,2) {mustBeReal, mustBeNonnegative, mustBeFinite} = [0,0] 
            end
                        
            [~,~,~,~,~] = this.transform();
            
            if intervallo == [0,0]
                [~,I] = max(this.amps);
                intervallo_up = I+40;
                intervallo_do = I-40;         
            else
                intervallo = sort(intervallo);
                intervallo_up = intervallo(1);
                intervallo_do = intervallo(2);     
            end
                
            % Seleziona solo dati nell'intorno del picco selezionato 
            tmp_freq = this.frequencies(intervallo_do:intervallo_up);
            tmp_amps = this.amps(intervallo_do:intervallo_up);

            % Calcolo media e sigma       
            densita_probabilita_freq = tmp_amps/sum(tmp_amps);
            frequenze_media = sum(tmp_freq .* densita_probabilita_freq);
            frequenzequadre_media =  sum((tmp_freq.^2).* densita_probabilita_freq);
            frequenze_var = frequenzequadre_media - frequenze_media^2;
            frequenze_sigma = sqrt(frequenze_var);
        end        
    end

    methods (Hidden)
        function [fig, ax] = plotTransform(this, yData, fileName, isAbsPlot)
            arguments
                this,               
                yData (:,1) double {mustBeReal, mustBeFinite},
                fileName (1,1) string = "",
                isAbsPlot (1,1) logical = 1,
            end
                                               
            fig = figure('units','inch','position',[0,0,this.figureWidth,this.figureHeight]);
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
            
            % Esporta figura
            if (strlength(fileName) > 0)
                exportFigure(fig, ax, fileName,this.fontSize, this.figureWidth, this.figureHeight);
            end                   
        end
    end
end
