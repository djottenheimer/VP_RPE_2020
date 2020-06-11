function loadedData = loadData_ott(task)

master_root = ottBari2020_root();
switch task
    case 'intBlocks'
        load(fullfile(master_root, 'Data', 'Modeling', 'DataForModeling', 'ModData_intBlocks.mat'), 'CS')
        loadedData = CS;
    case 'threeOutcomes'
        load(fullfile(master_root, 'Data', 'Modeling', 'DataForModeling', 'ModData_TH.mat'), 'TH')
        loadedData = TH;
    case 'cue'
        load(fullfile(master_root, 'Data', 'Modeling', 'DataForModeling', 'ModData_Cued.mat'), 'CD')
        loadedData = CD;
    otherwise
        error('task not found\n')
end