function [errRMS,errVar,meanShift] = fcn_Points_calcPairStatistics(pairedXYdata)
% fcn_Points_calcPairStatistics
% Calculates the statistics of the errors between two sets of paired XY
% data points and an overall X,Y map shift to match the centroids of the
% data sets.
%
%
% FORMAT:
%
%       [errRMS,errVar,meanShift = fcn_Points_calcPairStatistics(pairedXYdata)
%
% INPUTS:
%
%      pairedXYdata: an N x 4 matrix with [X1 Y1 X2 Y2] data in each row.
%
% OUTPUTS:
%
%      errRMS: the RMS distances associated with the point-to-point match errors
%      errVar: the variance in the point-to-point match error distances
%      meanShift: the [X Y] shift that would cause the centroids of the
%                 two data sets to match
%
% DEPENDENCIES
%
%      none
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_Points_calcPairStatisticss.m for a full
%       test suite.
%
% This function was written on 2022_01_21 by C. Beal
% Questions or comments? cbeal@bucknell.edu

% Revision history:
%     2022_01_21
%       -- wrote the code

% Compute a vector of distances between paired points
errDist = sqrt((pairedXYdata(:,3)-pairedXYdata(:,1)).^2 + (pairedXYdata(:,4) - pairedXYdata(:,2)).^2);
% Compute the RMS error distance between paired points
errRMS = sqrt(mean(errDist.^2));
% Compute the variance in the error distance between paired points
errVar = std(errDist)^2;

% Compute a row vector with the mean x,y shift between the centroids of
% the two data sets
meanShift = [mean(pairedXYdata(:,3)-pairedXYdata(:,1),1), mean(pairedXYdata(:,4)-pairedXYdata(:,2),1)];

end

