classdef fourierTransform < handle
    properties
        data (:,1) double {mustBeReal, mustBeFinite}
        dt (1,1) double {mustBeReal, mustBeFinite}
        frequencies (:,1) double {mustBeReal, mustBeFinite}
        amps (:,1) double {mustBeReal, mustBeFinite}
        phases (:,1) double {mustBeReal, mustBeFinite}
        dF (1,1) double {mustBeReal, mustBeFinite}
        name (1,1) string
        color (1,3) double {mustBeReal, mustBeFinite}
        freqlim (1,2) double {mustBeReal, mustBeFinite}
        mirrored (1,1) logical
        omegasAsXAxis (1,1) logical 
        showDots (1,1) logical 
        showLines (1,1) logical 
    end


    methods
        
        % -----------------------------------------------------------------
        
        function this = fourierTransform()
            this.data = [];
            this.dt = 1;
            this.frequencies = [];
            this.dF = 0;
            this.name = "FFT";
            this.color = [0 0 1];
            this.mirrored = 0;
            this.freqlim = [0,0];
            this.omegasAsXAxis = 0;
            this.showDots = 1;
            this.showLines = 1;
        end

        % -----------------------------------------------------------------
        
        function [frequencies, amps, phases, dF] = transform(this)
            arguments
                this                
            end
            
            % frequenza massima                      
            Fs = 1/this.dt;                                
            % Incertezza sulle frequenze
            dF = Fs/length(this.data); 
            % Vettore frequenze e pulsazioni
            this.frequencies = (-Fs/2:dF:Fs/2-dF)';
            % FFT dati
            fft_data = fftshift(fft(this.data));
            
            % Vettore Ampiezze e Vettore fasi 
            this.amps = abs(fft_data)/length(this.data);
            this.phases = angle(fft_data) / pi;
            % Output
            amps = this.amps;
            phases = this.phases;
            frequencies = this.frequencies;
        end
        
        % -----------------------------------------------------------------
        
        function [frequencies, phases, amps, dF] = plotTransform(this, fileName)
            [frequencies, phases, amps, dF] = plotAbsTransform(this, fileName);
        end

        % -----------------------------------------------------------------
        
        function [frequencies, phases, amps, dF] = plotAbsTransform(this, fileName)
            arguments
                this,
                fileName (1,1) string = ""
            end
            
            % Trasformata
            [frequencies, amps, phases, dF] = this.transform();
                        
            fig = figure();
            ax = axes();
            hold on
            
            % Scegli se mostrare frequenze o pulsazioni
            if this.omegasAsXAxis
                xAxisData = frequencies*2*pi;
                dX = dF*2*pi;
                xLimits = this.freqlim*2*pi;
                xlabel("Pulsazione [Hz]");
            else
                xAxisData = frequencies;
                dX = dF;
                xLimits = this.freqlim;
                xlabel("Frequenza [Hz]");
            end

            % Disegna punti righe verticali
            if this.showDots
                scatter(xAxisData,amps,100,"blue",".");
            end
            % Disegna righe verticali
            if this.showLines
                for ii = 1:length(frequencies)
                    line([xAxisData(ii) xAxisData(ii)],[0 amps(ii)],"Color","blue");
                end
            end
            
            % Imposta limiti
            if this.freqlim == [0,0] 
                if this.mirrored
                    xlim([min(xAxisData) - dX, max(xAxisData) + dX]);
                else
                    xlim([0, max(xAxisData) + dX]);
                end
            else
                xlim(sort(xLimits));
            end
                       
            ylabel("Modulo [V]");
            title(this.name);           
            ax.FontSize = 14;
            grid minor    

            if (strlength(fileName) > 0)
                exportFigure(fig, ax, fileName);
            end
                    
        end

        % -----------------------------------------------------------------

        % -----------------------------------------------------------------
        
        function [frequencies, phases, amps, dF] = plotPhaseTransform(this, fileName)
            arguments
                this,
                fileName (1,1) string = ""
            end
            
            % Trasformata
            [frequencies, amps, phases, dF] = this.transform();
                        
            fig = figure();
            ax = axes();
            hold on
            
            % Scegli se mostrare frequenze o pulsazioni
            if this.omegasAsXAxis
                xAxisData = frequencies*2*pi;
                dX = dF*2*pi;
                xLimits = this.freqlim*2*pi;
                xlabel("Pulsazione [Hz]");
            else
                xAxisData = frequencies;
                dX = dF;
                xLimits = this.freqlim;
                xlabel("Frequenza [Hz]");
            end

            % Disegna punti righe verticali
            if this.showDots
                scatter(xAxisData,phases,100,"blue",".");
            end
            % Disegna righe verticali
            if this.showLines
                for ii = 1:length(frequencies)
                    line([xAxisData(ii) xAxisData(ii)],[0 phases(ii)],"Color","blue");
                end
            end
            
            % Imposta limiti
            if this.freqlim == [0,0] 
                if this.mirrored
                    xlim([min(xAxisData) - dX, max(xAxisData) + dX]);
                else
                    xlim([0, max(xAxisData) + dX]);
                end
            else
                xlim(sort(xLimits));
            end

            ylim([-1.1 1.1]);
                       
            ylabel("Fase/\pi");
            title(this.name);           
            ax.FontSize = 14;
            grid minor    

            if (strlength(fileName) > 0)
                exportFigure(fig, ax, fileName);
            end
                    
        end

        % -----------------------------------------------------------------

    end
end
