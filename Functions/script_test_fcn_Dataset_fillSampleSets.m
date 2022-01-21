% script_test_fcn_Dataset_fillSampleSets.m
% This is a script to exercise the function:
% fcn_Path_fillSampleSets.m

% This function was adapted on 2022_01_12 by C. Beal from S. Brennan's
% script_test_fcn_Path_fillSamplePaths
% Questions or comments? sbrennan@psu.edu

% Revision history
%     2022_01_12
%     -- adapted for XY datasets without path ordering

close all
clc

% Call the function to fill in an array of "dataset" type
datasets_array = fcn_Dataset_fillSampleSets;

% We can even save one of these as a single "dataset"
single_path = datasets_array{1};

% % Convert paths to traversals structures. Each traversal instance is a
% % "traversal" type, and the array called "data" below is a "traversals"
% % type.
% for i_Path = 1:length(paths_array)
%     traversal = fcn_Path_convertPathToTraversalStructure(paths_array{i_Path});
%     data.traversal{i_Path} = traversal;
% end
%

% Call the plot command to show results in XY
% Show result
figure
hold on
for i_set = 1:length(datasets_array)
    plot(datasets_array{i_set}(:,1),datasets_array{i_set}(:,2),'.','Markersize',20);
    text(datasets_array{i_set}(1,1),datasets_array{i_set}(1,2),'Start');
end


