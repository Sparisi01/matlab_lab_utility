function ann = textBox(text, placement, ax, fontsize, hor_padding, ver_padding)
    % textBox(text, placement, ax, fontsize, hor_padding, ver_padding)
    %
    % text          string
    % placement     string = "northeast"
    % ax                   = gca
    % fontsize      double = 14
    % hor_padding   double = 0.015
    % ver_padding   double = 0.015
    %
    % All arguments are optional except for 'text'. The 'placement' argument take one 
    % of the following values: "northeast", "northwest", "southeast" or "southwest".
    %
    % A minimal working example is: 
    %   textBox("Hello, World!");

    arguments
        text string
        placement string = "northeast"
        ax = gca
        fontsize double {mustBeFinite, mustBeReal, mustBePositive} = 14
        hor_padding double {mustBeFinite, mustBeReal} = 0.015
        ver_padding double {mustBeFinite, mustBeReal} = 0.015
    end

    axes(ax);
    ann = annotation("textbox", [0,0,0.2,0.2], ...
        "BackgroundColor", [1,1,1], ...
        "FontSize", fontsize, ...
        "String", text, ...
        'FitBoxToText', 'on' ...
    );
    pause(0);

    ann = geograficPlacement(ann, placement, ax, hor_padding, ver_padding);
end