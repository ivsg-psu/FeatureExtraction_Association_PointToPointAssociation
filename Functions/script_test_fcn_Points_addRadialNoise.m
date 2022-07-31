% script_test_fcn_points_addRadialNoise
% This is a script to exercise the function:
% fcn_points_addRadialNoise.m

% This script was written on 2022_07_18 by Shashank Garikipati
% Questions or comments?

% Revision history
%     2022_07_18
%     -- wrote the code

% Clear out the workspace
close all
clearvars

% Load some test data sets
origXYdatasets = fcn_Points_fillPointSampleSets;

%% Example 1
% shows addition of radial noise

new_set1 = {[0.0 0.0]};
for i =  1:100
    datasets_out = fcn_Points_addRadialNoise(new_set1,0.5);
    fcn_Points_plotSetsXY(new_set1,1);
    fcn_Points_plotSetsXY(datasets_out,1);
    axis([-1 1 -1 1])
    axis square
end

%% Example 2
% Addition of radial noise for one data set

new_set2 = {origXYdatasets{1}}; %#ok<CCAT1> 
datasets_out2 = fcn_Points_addRadialNoise(new_set2,1);
fcn_Points_plotSetsXY(new_set2,2);
fcn_Points_plotSetsXY(datasets_out2,2);


%% Example 3
% Addition of radial noise for one data set

new_set3 = origXYdatasets;
datasets_out3 = fcn_Points_addRadialNoise(new_set3,1);
fcn_Points_plotSetsXY(new_set3,2);
fcn_Points_plotSetsXY(datasets_out3,2);

%% Example 4
% association of cone after adding radial noise

cone_location = {[1,1]};
datasets_out4 = fcn_Points_addRadialNoise(cone_location,1);
[pairedXYdata, numMatches, nonMatchesA]  = fcn_Points_pairXYdata(cone_location{1},datasets_out4{1},1);

figure(4);
clf
hold on
grid on 
for i = 1:numMatches
    plot([pairedXYdata(i,1) pairedXYdata(i,3)],[pairedXYdata(i,2) pairedXYdata(i,4)],'o','MarkerSize',10,'MarkerFaceColor',[1 .647 .25])
end
quiver(pairedXYdata(:,1),pairedXYdata(:,2),pairedXYdata(:,3)-pairedXYdata(:,1),pairedXYdata(:,4)-pairedXYdata(:,2),0,'k','linewidth',1)
% Also plot the non-matches in data set A with circles in various colors
for i = 1+numMatches:nonMatchesA+numMatches
    plot(pairedXYdata(i,1),pairedXYdata(i,2),'o','MarkerSize',10,'MarkerFaceColor',[1 .647 .25])
end
% Lastly, plot the non-matches in data set B with squares in various colors
for i = 1+numMatches+nonMatchesA:nonMatchesA+numMatches+nonMatchesA
    plot(pairedXYdata(i,3),pairedXYdata(i,4),'o','MarkerSize',10,'MarkerFaceColor',[1 .647 .25])
end

%% Example 5
% association of cone after adding radial noise, the cone may or may not be
% paired depending on its location. random noise of upto 2 units, but
% max association of only 1 unit

cone_location = {[1,1]};
datasets_out4 = fcn_Points_addRadialNoise(cone_location,2);
[pairedXYdata, numMatches, nonMatchesA,]  = fcn_Points_pairXYdata(cone_location{1},datasets_out4{1},1);

figure(5);
clf
hold on
grid on
for i = 1:numMatches
    plot([pairedXYdata(i,1) pairedXYdata(i,3)],[pairedXYdata(i,2) pairedXYdata(i,4)],'o','MarkerSize',10,'MarkerFaceColor',[0 1 0])
end
quiver(pairedXYdata(:,1),pairedXYdata(:,2),pairedXYdata(:,3)-pairedXYdata(:,1),pairedXYdata(:,4)-pairedXYdata(:,2),0,'k','linewidth',1)
% Also plot the non-matches in data set A with circles in various colors
for i = 1+numMatches:nonMatchesA+numMatches
    plot(pairedXYdata(i,1),pairedXYdata(i,2),'o','MarkerSize',10,'MarkerFaceColor',[0 0 1])
end
% Lastly, plot the non-matches in data set B with squares in various colors
for i = 1+numMatches+nonMatchesA:nonMatchesA+numMatches+nonMatchesA
    plot(pairedXYdata(i,3),pairedXYdata(i,4),'o','MarkerSize',10,'MarkerFaceColor',[1 0 0])
end

%% Example 6
% association of cone after adding radial noise, the cone may or may not be
% paired depending on its location. random noise of upto 2 units, but
% max association of only 1 unit

cone_location = {[1,1;0,0;1,2; 4,3; 5,7;3,3;3,5; 5, 4;1,5]};
datasets_out4 = fcn_Points_addRadialNoise(cone_location,1.5);

[pairedXYdata, numMatches, nonMatchesA, nonMatchesB]  = fcn_Points_pairXYdata(cone_location{1},datasets_out4{1},1);

figure(6);
clf
hold on
grid on

for i = 1:numMatches
    plot([pairedXYdata(i,1) pairedXYdata(i,3)],[pairedXYdata(i,2) pairedXYdata(i,4)],'o','MarkerSize',10,'MarkerFaceColor',[0 1 0])
end
quiver(pairedXYdata(:,1),pairedXYdata(:,2),pairedXYdata(:,3)-pairedXYdata(:,1),pairedXYdata(:,4)-pairedXYdata(:,2),0,'k','linewidth',1)
% Also plot the non-matches in data set A with circles in various colors
for i = 1+numMatches:nonMatchesA+numMatches
    plot(pairedXYdata(i,1),pairedXYdata(i,2),'o','MarkerSize',10,'MarkerFaceColor',[1 0 0])
end
% Lastly, plot the non-matches in data set B with squares in various colors
for i = 1+numMatches+nonMatchesA:nonMatchesA+numMatches+nonMatchesA
    plot(pairedXYdata(i,3),pairedXYdata(i,4),'o','MarkerSize',10,'MarkerFaceColor',[0 0 1])
end

