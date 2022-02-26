function [condMinVal,condMinInd] = conditionalMin(valueVect,conditionalVect,conditional)
% conditionalMin
% Finds the minimum value element in a vector of values where the element
% in an associated conditional vector is true. The value vector and
% conditional vector *may* be the same if the conditional is based on the
% values in the value vector (finding the minimum of non-negative values,
% for example).
%
% INPUTS:
%
%      valueVect: A vector containing values over which to find the min.
%      conditionalVect: A vector containing values or flags against which
%           the conditional expression will be tested.
%      conditional: A string which assumes that the conditional vector will
%           be prepended as the first word and defines a conditional test
%           (e.g. input '== -1' expands to 'conditionalVect == -1'). This
%           ensures that vectors of any name can be passed into the
%           function and used for the conditional test.
%
% OUTPUTS:
%
%      minVal: the minimum value of the valueVector that satisfies the
%           condition
%      minInd: the position of the minimum value in valueVector that
%           satisfies the condition
%
% DEPENDENCIES:
%
%      none
%
% EXAMPLES:
%
%       [minVal,minInd] = conditionalMin(valueVector,flagVector,'== 0')
%        - Finds the minimum value of valueVector where the associated
%          flag in flagVector is zero. Elements in valueVector where the
%          associated element in flagVector is non-zero are ignored.
%
% This function was written by C. Beal
% Questions or comments? cbeal@bucknell.edu

% Revision history:
%     2022_02_26
%     -- wrote the code 

% Check for size agreement
if size(valueVect) ~= size(conditionalVect)
    % If vectors are the same size but oriented differently, fix and
    % continue but warn the user
    if size(valueVect,1) > size(conditionalVect,1) && size(valueVect,1) == size(conditionalVect,2)
        warning('Input vectors are column and row vectors, respectively. Reshaping and continuing.')
        conditionalVect = conditionalVect';
    elseif size(valueVect,1) < size(conditionalVect,1) && size(valueVect,1) == size(conditionalVect,2)
        warning('Input vectors are row and column vectors, respectively. Reshaping and continuing.')
        minVect = minVect';
    % If the vectors are truly different sizes, error and stop
    else
        error('Incompatible input vector sizes')
    end
end

% Find the indices of elements of the vector that satisfy the condition
condInds = find(eval(['conditionalVect' conditional]));

% Find the minimum of the sub-vector of the elements that satisfy the condition
[minVal,minInd] = min(valueVect(condInds));

% The conditional minimum is the same as the minimum, but the index must be
% adjusted to address the element in the original vector, not the
% sub-vector with valid elements
condMinVal = minVal;
condMinInd = condInds(minInd);

