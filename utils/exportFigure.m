function exportFigure(image_figure, image_axes, image_name, fontSize, image_width, image_height)
    % -------------------------------------------------
    % Funzione per l'export automatico di una figura
    % in formato png. 
    % -------------------------------------------------

    arguments
        image_figure,
        image_axes,
        image_name,
        fontSize (1,1) double {mustBeReal, mustBeFinite} = 14
        image_width (1,1) double {mustBeReal, mustBeFinite} = 8; % 8 inches
        image_height (1,1) double {mustBeReal, mustBeFinite} = 6; % 6 inches -> ratio 4:3
    end

    for ax = image_axes 
        set(ax, "FontSize", fontSize);
    end

    set(image_figure, 'PaperUnits', 'inches');
    set(image_figure, 'PaperSize', [image_width image_height]);

    set(image_figure,'InvertHardcopy','on');
    set(image_figure,'PaperUnits', 'inches');
    set(image_figure,'PaperPosition', [0, 0, image_width, image_height]);

    % Save the file as PNG
    print(image_name,'-dpng','-r300');

end