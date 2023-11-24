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

    NORTH = 1; SOUTH = 2; WEST = 3; EAST = 4;
    vertical_alignement = NORTH;
    horizontal_alignement = EAST;
    if placement == "northeast"
        vertical_alignement = NORTH;
        horizontal_alignement = EAST;
    elseif placement == "northwest"
        vertical_alignement = NORTH;
        horizontal_alignement = WEST;
    elseif placement == "southeast"
        vertical_alignement = SOUTH;
        horizontal_alignement = EAST;
    elseif placement == "southwest"
        vertical_alignement = SOUTH;
        horizontal_alignement = WEST;
    else
        error("The argument 'placement' should have one of the following values: 'northeast', 'northwest', 'southeast', 'southwest'. Instead '" + placement + "' was given.");
    end

    x = 1; y = 2; w = 3; h = 4;
    if horizontal_alignement == EAST
        ann.Position(x) = ax.Position(x) + ax.Position(w) - ann.Position(w) - hor_padding;
    elseif horizontal_alignement == WEST
        ann.Position(x) = ax.Position(x) + hor_padding;
    end
    if vertical_alignement == NORTH
        ann.Position(y) = ax.Position(y) + ax.Position(h) - ann.Position(h) - ver_padding;
    elseif vertical_alignement == SOUTH
        ann.Position(y) = ax.Position(y) + ver_padding;
    end
end