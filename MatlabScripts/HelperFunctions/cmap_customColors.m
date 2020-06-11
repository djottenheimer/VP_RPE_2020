function c = cmap_customColors(m, color)    

switch color
    case 'vermillion'
        c(:, 1) = linspace(0.8, 1, m);
        c(:, 2) = linspace(0.4745, 1, m);
        c(:, 3) = linspace(0.6549, 1, m);
    case 'skyBlue'
        c(:, 1) = linspace(0.3373, 1, m);
        c(:, 2) = linspace(0.7059, 1, m);
        c(:, 3) = linspace(0.9137, 1, m);
    case 'whiteBlue'
        c(:, 1) = linspace(0.9, 0, m);
        c(:, 2) = linspace(0.9447, 0.4470, m);
        c(:, 3) = linspace(0.9741, 0.7410, m);
    case 'grays'
        c(:, 1) = linspace(0.1, 0.9, m);
        c(:, 2) = linspace(0.1, 0.9, m);
        c(:, 3) = linspace(0.1, 0.9, m);
    otherwise
        error('Choose appropriate color name')
end