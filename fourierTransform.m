classdef fourierTransform < handle
    properties
        data (:,1) double {mustBeReal, mustBeFinite}
        frequencies (:,1) double {mustBeReal, mustBeFinite}
        dt (1,1) double {mustBeReal, mustBeFinite}        
        amps (:,1) double {mustBeReal, mustBeFinite}
        phases (:,1) double {mustBeReal, mustBeFinite}
        dF (1,1) double {mustBeReal, mustBeFinite}        
        name (1,1) string
        xLabel (1,1) string
        yLabel_phase (1,1) string
        yLabel_abs (1,1) string
        color (1,3) double {mustBeReal, mustBeFinite}
        xAxisLim (1,2) double {mustBeReal, mustBeFinite}
        xAxisAsOmegas (1,1) logical       
        mirrored (1,1) logical
        tollerance (1,1) double {mustBeReal, mustBeFinite}
    end


    methods
        
        % -----------------------------------------------------------------
        
        function this = fourierTransform()
            this.data = [];
            this.frequencies = [];
            this.dt = 1;
            this.dF = 0;
            this.name = "FFT";
            this.xLabel = "[Hz]";
            this.yLabel_phase = "\angle{FFT} [deg]";
            this.yLabel_abs = "|FFT| [V]";
            this.color = [0 0 1];
            this.mirrored = 0;
            this.xAxisLim = [0,0];
            this.xAxisAsOmegas = 0;
            this.tollerance = 1e-5;
        end

        % -----------------------------------------------------------------
        % Utilizzo della funzione FFT per estrarre vettore di ampiezze e fasi dal vettore dati
        % https://www.gaussianwaves.com/2015/11/interpreting-fft-results-obtaining-magnitude-and-phase-information/           
        function [frequencies, amps, phases, dF] = transform(this)
            arguments
                this                
            end
                     
            % FFT dati
            fft_data = fftshift(fft(this.data));
            % Frequenza massima risolvibile                      
            Fs = 1/this.dt;                                
            % Incertezza sulle frequenze
            this.dF = Fs/length(this.data); 
           
            % Vettore frequenze
            this.frequencies = (-Fs/2:this.dF:Fs/2-this.dF)';
                       
            % Vettore Ampiezze
            this.amps = abs(fft_data)/length(this.data);
            
            % Vettore Fasi   
            tmp_fft_data = fft_data;   
            threshold = max(abs(fft_data))*this.tollerance;
            tmp_fft_data(abs(fft_data)<threshold) = 0;  
            this.phases = atan2(imag(tmp_fft_data),real(tmp_fft_data))*180/pi; 
                
            % Output
            amps = this.amps;
            phases = this.phases;
            frequencies = this.frequencies;
            dF = this.dF;
           
        end
        
        % -----------------------------------------------------------------
        
        function [frequencies, phases, amps, dF, fig, ax] = plotAbsTransform(this, fileName)
            arguments
                this,
                fileName (1,1) string = ""
            end
                        
            [frequencies, amps, phases, dF] = this.transform();

            [fig, ax] = plotTransform(this, amps, fileName, 1);
        end

        % -----------------------------------------------------------------
        
        function [frequencies, phases, amps, dF , fig, ax] = plotPhaseTransform(this, fileName)
            arguments
                this,
                fileName (1,1) string = ""
            end
                        
            [frequencies, amps, phases, dF] = this.transform();
            [fig, ax] = plotTransform(this, phases, fileName, 0);
        end

        % -----------------------------------------------------------------
        
    end
    methods (Hidden)
        function [fig, ax] = plotTransform(this, yData, fileName, isAbsPlot)
            arguments
                this,               
                yData (:,1) double {mustBeReal, mustBeFinite},
                fileName (1,1) string = "",
                isAbsPlot (1,1) logical = 1,
            end
                                               
            fig = figure();
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

            if (strlength(fileName) > 0)
                exportFigure(fig, ax, fileName);
            end
                    
        end

        % -----------------------------------------------------------------
    end
end
