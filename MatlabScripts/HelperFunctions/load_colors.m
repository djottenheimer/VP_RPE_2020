function Colors = load_colors()
colors = {'sucrose',[1  0.65  0.1];
    'maltodextrin',[.9  0.3  .9];
    'water',[0.00  0.75  0.75];
    'total',[0.3  0.1  0.8];
    'exc',[0 113/255 188/255];
    'inh',[240/255 0 50/255];
    'NAc',[0.5  0.1  0.8];
    'VP',[0.3  0.7  0.1];
    'block1',[143  158  173]/255; %1922 in color book
    'block2',[204  50  69]/255;
    'block_total',[74  31  60]/255;
    'Cue',[38 56 137]/255; %1976 in color book
    'PE',[146 89 159]/255;
    'RD',[229 36 40]/255;
    'blockrose', [0.3  0.85  0.7];
    %'blockrose', [0.3  0.65  0.3];
    'blockrodextrin', [0.25  0.25  1];
    'gray', [.8 .8 .8];
    'sucrosesnp',[.95  0.55  0.15];
    'maltodextrinsnp',[.9  0.3  .9];
    'sucrosemnp',[.98  0.8  0.35];
    'maltodextrinmnp',[.9  0.75  1];
    'sucrosesboth',[.8  0.7  0.15];
    'maltodextrinsboth',[.7  0.45  .9];
    'sucrosemboth',[.95  0.85  0.15];
    'maltodextrinmboth',[0.85  0.65  1];
    'mean',[0.6 0.6 0.6];
    'current',[0 0.25 0.25];
    'rpe',[0 0.9 0.5];
    'cueeffect',[0.7 0.25 0.2];
    };

Colors = containers.Map(colors(:,1),colors(:,2));
end