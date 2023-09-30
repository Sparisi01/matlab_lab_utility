function exportFigure(image_figure, image_axes, image_name)
% image_figure = gcf
% image_axes = gca
% image_name = "percorso/nome_del_file"
%
% Esporta un'immagine di proporzione 4:3 con la dimensione del font di 14
%
% E' preferibile non esporre i parametri per le dimensioni in modo
% da fissare uno stile omogeneo tra i grafici.
for ax = image_axes 
    set(ax, "FontSize", 14);
end

image_width = 8; % 8 inches
image_height = 6; % 6 inches -> ratio 4:3
set(image_figure, 'PaperUnits', 'inches');
set(image_figure, 'PaperSize', [image_width image_height]);

set(image_figure,'InvertHardcopy','on');
set(image_figure,'PaperUnits', 'inches');
set(image_figure,'PaperPosition', [0, 0, image_width, image_height]);

% Save the file as PNG
print(image_name,'-dpng','-r300');

end