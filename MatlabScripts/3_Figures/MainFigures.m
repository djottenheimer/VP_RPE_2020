%% Paper Figures
% This script calls most others such that you can generate all main
% figures by running it (excluding optogenetics experiments)
%% Reset
clearvars;
close all;
load ('RAWintBlocks.mat');
load ('R_intBlocks.mat');
load ('RAWCued');
load ('R_cued_all');
load ('RAWTH');
load ('R_TH');
load ('RewHist2R');
load ('RewHistCued');
load ('RewHistTH');
load ('evolution_data.mat')
load ('evolution_info.mat')
Colors = load_colors();

%% Behavior
% Produces most of the behavior panels in the paper
disp('Behavior')
plot_sub_figures = true;
perform_stats = true;
plot_region = 'VP'; %'VP' 'NA' or 'Both'
BehavStats=all_behavior(RAWblocks, RAWCued, plot_region, plot_sub_figures, perform_stats, Colors);

%% Interspersed sucrose/maltodextrin

plot_value=1;
load_fits=1;
if load_fits
    load ('intBlocks_MLEfits.mat');
    bm_RD=select_RPEmods(os, 'RD', 'particularModels', {'mean','curr','base'},'plotmodels_Flag',false);
    bm_RD_BIC=select_RPEmods(os, 'RD', 'particularModels', {'mean','curr','base'},'plotmodels_Flag',false,'scoretouse','BIC');
    bm_cue=select_RPEmods(os, 'cue', 'particularModels', {'mean','curr','base'},'plotmodels_Flag',false);
end

disp('First experiment figures')
IntStats=plot_int(Colors, R_blocks, RewHist2R, bm_RD, RAWblocks, os, bm_cue, plot_value, bm_RD_BIC);

%% Three rewards (sucrose/maltodextrin/water)

load_fits=1;
if load_fits
    load ('threeOutcomes_MLEfits.mat');
    bm_RD=select_RPEmods(os, 'RD', 'particularModels', {'mean','curr','base'},'plotmodels_Flag',false);
    bm_RD_BIC=select_RPEmods(os, 'RD', 'particularModels', {'mean','curr','base'},'plotmodels_Flag',false,'scoretouse','BIC');
    bm_cue=select_RPEmods(os, 'cue', 'particularModels', {'mean','curr','base'},'plotmodels_Flag',false);
end


plot_sub_figures = true;
disp('Three outcomes figure')
ThreeStats=plot_TH(Colors, R_TH, RewHistTH, bm_RD, bm_cue, RAWTH, os, plot_sub_figures, bm_RD_BIC);

%% Blocks
%blocks data
disp('Blocks figure')
plot_sub_figures = true;

load_fits=1;
if load_fits
    load ('intBlocks_MLEfits.mat');
    bm_RD=select_RPEmods(os, 'RD', 'particularModels', {'mean','curr','base'},'plotmodels_Flag',false);
end

selectivity_time_periods = [0, 1;-2.6, 0; 0.4, 3]; %seconds, cue, PE, RD
R_blocks=plot_blocks(Colors, R_blocks, events, trial_info, plot_sub_figures ,bm_RD);

%% Cued figure
%summary of cue experiment
disp('Cued figure')

load_fits=1;
if load_fits
    load ('cue_MLEfits.mat');
    bm_RD=select_RPEmods(os, 'RD', 'particularModels', {'mean','mean_cue','curr','curr_cue','base','base_cue'},'plotmodels_Flag',false);
    bm_cue=select_RPEmods(os, 'cue', 'particularModels', {'mean','mean_cue','base','base_cue'},'plotmodels_Flag',false);    
end

ROI=[3 4 9 10];
DOI=11:20;
CuedStats=plot_cued(Colors, R_cued, RewHistCued, bm_RD, DOI, ROI, os, bm_cue, RAWCued);