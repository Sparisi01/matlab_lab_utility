function obj = geograficPlacement(obj, placement, ax, hor_padding, ver_padding)
    % geograficPlacement(obj, placement, ax, hor_padding, ver_padding)
    %
    % obj           class with 'Position' attribute
    % placement     string = "northeast"
    % ax                   = gca
    % hor_padding   double = 0.015
    % ver_padding   double = 0.015
    %
    % All arguments are optional except for 'obj'. The 'placement' argument take one 
    % of the following values: "northeast", "northwest", "southeast" or "southwest".

    arguments
        obj
        placement string = "northeast"
        ax = gca
        hor_padding double {mustBeFinite, mustBeReal} = 0.015
        ver_padding double {mustBeFinite, mustBeReal} = 0.015
    end

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
        obj.Position(x) = ax.Position(x) + ax.Position(w) - obj.Position(w) - hor_padding;
    elseif horizontal_alignement == WEST
        obj.Position(x) = ax.Position(x) + hor_padding;
    end
    if vertical_alignement == NORTH
        obj.Position(y) = ax.Position(y) + ax.Position(h) - obj.Position(h) - ver_padding;
    elseif vertical_alignement == SOUTH
        obj.Position(y) = ax.Position(y) + ver_padding;
    end
end