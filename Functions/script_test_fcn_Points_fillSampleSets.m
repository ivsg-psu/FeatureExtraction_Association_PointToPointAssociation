% script_test_fcn_Points_fillSampleSets.m
% This is a script to exercise the function:
% fcn_Points_fillSampleSets.m

% This function was adapted on 2022_01_12 by C. Beal from S. Brennan's
% script_test_fcn_Path_fillSamplePaths
% Questions or comments? sbrennan@psu.edu

%% Revision history
% 2022_01_12
% -- adapted for XY datasets without path ordering
% 2022_07_01
% -- added sections and two test examples

clc, close all
%% Setup

% Call the function to fill in an array of "dataset" type
datasets_array = fcn_Points_fillPointSampleSets;


%% Example 1
% We can save one of these as a single "dataset"

single_path = {datasets_array{1}};

fcn_Points_plotSetsXY(single_path,1);
text(single_path{1}(1,1),single_path{1}(1,2),'Start');


%% Example 2
% illustrating the different XY points on different plots

for i_set = 1:length(datasets_array)
    fcn_Points_plotSetsXY(datasets_array(i_set),(i_set+1));
    text(datasets_array{i_set}(1,1),datasets_array{i_set}(1,2),'Start');
end

%% Example 3
% All plots on the same graph

fcn_Points_plotSetsXY(datasets_array,5);


%%


