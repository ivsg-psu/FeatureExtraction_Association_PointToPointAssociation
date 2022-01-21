function [pairedXYdata, numMatches, nonMatchesA, nonMatchesB] = fcn_Dataset_pairXYdata(xyDataA,xyDataB,varargin)

% fcn_Dataset_pairXYdata
% Determines the closest pairs of matching X,Y data points from two
% different data sets of (potentially) different length. The algorithm uses
% a "mutual" pairing approach, where the nearest neighbor point is computed
% for each point in each data set and then the matches between the closest
% points are paired up. The matched pairs are listed in a single matrix,
% followed by unmatched pairs from data set A and then data set B. The
% unmatched data points are padded with NaN values. The number of matched
% points as well as number of unmatched points from data set A and data set
% B are returned.
%
%
% FORMAT:
%
%       [pairedXYdata, numMatches, nonMatchesA, nonMatchesB] = fcn_Dataset_pairXYdata(xyDataA,xyDataB,(maxDist))
%
% INPUTS:
%
%      xyDataA: an N x 2 matrix with [X Y] data in each row.
%      xyDataB: an N x 2 matrix with [X Y] data in each row.
%
%      (OPTIONAL INPUTS)
%      maxDist: a maximum radius outside of which to reject matched pairs
%
% OUTPUTS:
%
%      pairedXYdata: an N x 4 matrix with [X1 Y1 X2 Y2] data in each row
%      numMatches: the number of mutual matches in data sets A and B
%      nonMatchesA: the number of points in data set A without a mutual match
%      nonMatchesB: the number of points in data set B without a mutual match
%
% DEPENDENCIES
%
%      none
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_Dataset_pairXYdata.m for a full
%       test suite.
%
% This function was written on 2022_01_21 by C. Beal
% Questions or comments? cbeal@bucknell.edu

flag_do_debug = 0; % Flag to plot the results for debugging
flag_do_plots = 0; % Flag to plot the final results
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end

if flag_check_inputs == 1
    % Are there the right number of inputs?
    if nargin < 2 || nargin > 3
        error('Incorrect number of input arguments')
    end
    
    % Check the Dataset variables (not yet implemented fully)
    % fcn_Dataset_checkInputsToFunctions(Dataset, 'dataset');

end

% Does user want to show the plots?
if 3 == nargin
    maxDist = varargin{1};
else
    maxDist = inf;
end

% Make sure data sets are properly set up
if(size(xyDataA,2) ~= 2)
    error('XY Dataset A does not contain two coordinates in a column-oriented matrix');
elseif(size(xyDataB,2) ~= 2)
    error('XY Dataset B does not contain two coordinates in a column-oriented matrix');
end

% Binary compatibility analysis (mutually nearest neighbor points)
% First, associate data points in the sets using a nearest neighbor search
[idxAtoB,distAtoB] = knnsearch(xyDataB,xyDataA);
[idxBtoA,distBtoA] = knnsearch(xyDataA,xyDataB);

% Find mutual matches between nearest neighbors. These points are "most" in
% agreement in that they are mutual nearest neighbor matches. Also check
% the distance between matches to make sure they are appropriately matched
mutualA = zeros(max(length(idxAtoB),length(idxBtoA)),1);
for i = 1:length(idxAtoB)
    if(idxBtoA(idxAtoB(i)) == i && distAtoB(i) <= maxDist)
        mutualA(i) = 1;
    end
end

% Determine the number of matches made for output
numMatches = sum(mutualA);

% Determine the mutual match indices on the vector idxBtoA
mutualB = zeros(length(idxBtoA),1);
mutualB(idxAtoB(mutualA==1)) = 1;

% Organize the matches in a 4-column output matrix
pairedXYdata = [xyDataA(idxBtoA(idxAtoB(mutualA==1)),:), xyDataB(idxAtoB(mutualA==1),:)];

% Determine the number of non-matches in A
nonMatchesA = size(xyDataA,1) - numMatches;
nonMatchesB = size(xyDataB,1) - numMatches;

% Append the non-matched points at the end of the matrix
pairedXYdata = [pairedXYdata; xyDataA(mutualA == 0,:) nan(nonMatchesA,2)];
pairedXYdata = [pairedXYdata; nan(nonMatchesB,2) xyDataB(mutualB == 0,:)];

end

