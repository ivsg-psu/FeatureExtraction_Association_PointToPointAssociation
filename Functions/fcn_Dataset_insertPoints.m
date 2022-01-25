function patchStruct = fcn_Dataset_insertPoints(patchStruct,pointArray)
% fcn_Dataset_insertPoints
% Inserts a set of points into a patch structure, preserving the CCCW 
% ordering of the points such that a patch graphics object can be plotted
%
% FORMAT: 
%
%       updatedPatch = fcn_Dataset_plotPatch(patchStruct,pointArray)
%
% INPUTS:
%
%      patchStruct: a structure containing subfields of X and Y coordinates
%      in the following form
%           patchArray{i_patch}.X
%           patchArray{i_patch}.Y
%      Note that i_patch denotes an array of patch structures. Each 
%      structure will be plotted separately.
%      pointArray: an N x 2 matrix with column-wise X and Y coordinates 
%                  of points to insert into the patch
%
% OUTPUTS:
%
%      patchStruct: the updated patch structure with the points added
%
% DEPENDENCIES:
%
%      fcn_Dataset_determineAABB
%
%      ## NOT CURRENTLY USED: fcn_DataSet_checkInputsToFunctions
%
% EXAMPLES:
%      
%       See the script: script_test_fcn_Dataset_insertPoints.m for a full test
%       suite. 
%
% This function was written by C. Beal
% Questions or comments? cbeal@bucknell.edu 

% Revision history:
%     2022_01_25
%     -- wrote the code

flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end


%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _       
%  |_   _|                 | |      
%    | |  _ __  _ __  _   _| |_ ___ 
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |                  
%              |_| 
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_check_inputs == 1
    % Are there the right number of inputs?
    if nargin ~= 2
        error('Incorrect number of input arguments')
    end
    
    % Check to see that there were any points provided
    if isempty(pointArray)
        warning('Empty array of points to add, nothing to do.');
        return
    end
    
    % Check the point array input
    if size(pointArray,2) ~= 2
        error('Wrong format for point array input')
    end
end



%% Main body of the code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumPoints = size(pointArray,1);


% The first point will be the closest point, the highest index will be the
% one furthest around counter-clockwise, and the lowest index will be the
% most clockwise point. We want to find the most counter-clockwise
% half-space in which the test point fits.

% define the dot product for use in testing half-space membership
dot = @(x,y) x(:)'*y(:);

for i = 1:NumPoints
    % Find the closest patch point to the test point
    [closePt,~] = knnsearch([patchStruct.pointsX patchStruct.pointsY],pointArray(i,:),'nsmethod','exhaustive');
    prevPt = mod(closePt-2,length(patchStruct.pointsX))+1;
    nextPt = mod(closePt,length(patchStruct.pointsX))+1;
    dotProd = dot(pointArray(i,:) - [patchStruct.pointsX(closePt) patchStruct.pointsY(closePt)],[patchStruct.pointsX(nextPt) patchStruct.pointsY(nextPt)] - [patchStruct.pointsX(closePt) patchStruct.pointsY(closePt)]);
    if(0 < dotProd)
        if(1 == nextPt)
            patchStruct.pointsX = [patchStruct.pointsX(1:end); pointArray(i,1)];
            patchStruct.pointsY = [patchStruct.pointsY(1:end); pointArray(i,2)];            
        else
            patchStruct.pointsX = [patchStruct.pointsX(1:closePt); pointArray(i,1); patchStruct.pointsX(closePt+1:end)];
            patchStruct.pointsY = [patchStruct.pointsY(1:closePt); pointArray(i,2); patchStruct.pointsY(closePt+1:end)];
        end
    else
        if(1 == closePt)
            patchStruct.pointsX = [pointArray(i,1); patchStruct.pointsX(closePt:end)];
            patchStruct.pointsY = [pointArray(i,2); patchStruct.pointsY(closePt:end)];            
        else
            patchStruct.pointsX = [patchStruct.pointsX(1:prevPt); pointArray(i,1); patchStruct.pointsX(closePt:end)];
            patchStruct.pointsY = [patchStruct.pointsY(1:prevPt); pointArray(i,2); patchStruct.pointsY(closePt:end)];
        end
    end
end

% Update the axis-aligned bounding box (aabb)
patchStruct = fcn_Dataset_determineAABB(patchStruct);