function colorStruct = importColors_bb()
% Colors from Wong B. Nat Methods 8;6, 441 (2011).


workbookFile = fullfile(ottBari2020_root, 'MatlabScripts', 'HelperFunctions', 'colors.xlsx');

% Import the data
[~, ~, raw] = xlsread(workbookFile);
raw = raw(2:end, :);

charVectors = string(raw(:,[1, 2]));
charVectors(ismissing(charVectors)) = '';
charVectors = cellstr(charVectors);

for i = 1:size(charVectors, 1)
%     charVectors{i, 2} = str2double(charVectors{i, 2});
    
    color = str2num(charVectors{i, 2})./255;
    colorName = charVectors{i, 1};
    
    colorStruct.(colorName) = color;
end