% script_test_fcn_Points_plotSetsXY.m
% This is a script to exercise the function:
% fcn_Points_plotSetsXY.m
% This script was written on 2022_7_1 by Shashank Garikipati

%% Revision 
% -- Shashank Garikipati
% -- First write of code

clc, close all
%% Example 1
% Plots XY positions of one dataset
fig_num = 1;

dataset = {};

dataset{1} = [1,21;
              2,2;
              3,3;
              4,1;
              5,4;
              6,32;
              7,0.5];

fcn_Points_plotSetsXY(dataset,fig_num)

%% Example 2
% Plots XY positions of three datasets

dataset2 = {};
dataset2{1} = [1,21;
              2,2;
              3,3;
              4,1;
              5,4;
              6,32;
              7,0.5];

dataset2{2} = dataset{1}*1.8+0.3;

dataset2{3} = [9,2;
              8,0;
              2,5;
              1,1;
              6,4;
              10,31;
              2,1];

fcn_Points_plotSetsXY(dataset2,2)

%% Example 3
% Plots from a sample data set inside fcn_Points_fillPointSampleSets

fig_num = 3;
datasets3 = fcn_Points_fillPointSampleSets;

fcn_Points_plotSetsXY(datasets3,fig_num);

for i_set = 1:length(datasets3)
    text(datasets3{i_set}(1,1),datasets3{i_set}(1,2),'Start');
end

