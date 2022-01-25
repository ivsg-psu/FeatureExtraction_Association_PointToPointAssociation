function [patchArray, varargout] = fcn_Dataset_inferPrimitive(patchArray,varargin)
% fcn_Dataset_inferPrimitive
% Infer the most appropriate 2D primitive shape (rectangle, circle, or
% irregular). This may be extended into 3d in the future with truncated 
% rectangular pyramids and truncated cones as the primitives.
%
% FORMAT:
%
%       [patchArray, varargout] = fcn_Dataset_inferPrimitive(patchArray,{fig_num},{indices})
%
% INPUTS:
%
%      patchStruct: a structure containing subfields of X and Y coordinates
%      in the following form
%           patchArray{i_patch}.X
%           patchArray{i_patch}.Y
%      Note that i_patch denotes an array of patch structures. Each
%      structure will be plotted separately.
%
%
%      (OPTIONAL INPUTS)
%      fig_num: a figure number to plot into
%      indices: a vector of indices indicating which patch objects to plot
%
% OUTPUTS:
%
%      patchStruct: the updated patch structure with the AABB determined or
%      updated
%
%      (OPTIONAL OUTPUTS)
%      hpr: a vector of handles to the plotted primitives in order of the 
%           indices (if provided)
%
% DEPENDENCIES:
%
%      fcn_Dataset_determineAABB
%
%      ## NOT CURRENTLY USED: fcn_DataSet_checkInputsToFunctions
%
% EXAMPLES:
%
%       See the script: script_test_fcn_Dataset_inferPrimitive.m for a full test
%       suite.
%
% This function was written by C. Beal
% Questions or comments? cbeal@bucknell.edu

% Revision history:
%     2022_01_25
%     -- wrote the code

flag_do_debug = 0; % Flag to add debugging code
flag_check_inputs = 1; % Flag to perform input checking
flag_plot_prim = 0; % Flag to plot the axis-aligned bounding boxes

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
    if nargin < 1 || nargin > 3
        error('Incorrect number of input arguments')
    end
    
    % Check the data input
    
    % fcn_Path_checkInputsToFunctions(traversals, 'traversals');
    
end

% Did the user provide a figure number?
if 1 < nargin
    fig_num = varargin{1};
    figure(fig_num);
    flag_plot_prim = 1;
end

% Did the user provide a specific index vector?
if 3 == nargin
    idxVec = varargin{2};
else    
    idxVec = (1:length(patchArray))';
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

% Determine the axis-aligned bounding box (aabb) for each patch structure
NumPatches = length(idxVec);
% Loop through all of the patches
for i_patch = 1:NumPatches
    % Create an anonymous function to optimize for the circle center and
    % radius. This function is specific to each set of points for a patch
    centerFun = @(x) sum((x(3)^2 - ((patchArray(idxVec(i_patch)).pointsX-x(1)).^2 + (patchArray(idxVec(i_patch)).pointsY-x(2)).^2)).^2);
    % Make sure the axis-aligned bounding box is populated
    if isempty(patchArray(idxVec(i_patch)).aabb)
        % Update the axis-aligned bounding box (aabb)
        patchArray(idxVec(i_patch)) = fcn_Dataset_determineAABB(patchArray(idxVec(i_patch)));
    end
    % Set the initial guess for the center of the circle to the center of the bounding box
    x0 = [mean(patchArray(idxVec(i_patch)).aabb([1 3])); mean(patchArray(idxVec(i_patch)).aabb([2 4]))];
    % Set the initial guess for the radius of the circle to average of the
    % width and height of the bounding box
    x0(3) = mean([patchArray(idxVec(i_patch)).aabb(3)-patchArray(idxVec(i_patch)).aabb(1); patchArray(idxVec(i_patch)).aabb(4)-patchArray(idxVec(i_patch)).aabb(2)]);
    % Find the optimal center of the circle
    [soln,~] = fminsearch(centerFun,x0);
    % Compute the optimal radius for the circle
    patchArray(idxVec(i_patch)).primparams(1:2) = [soln(1) soln(2)];
    patchArray(idxVec(i_patch)).primparams(3) = soln(3);
end

% Find the best rectangle fit to each patch
% Assumption: 
% - Points are sorted such that they are adjacent to each other in a
%   counter-clockwise convention
% Algorithm:
% - Determine the farthest two data points from each other
% - Make these the extremum points of two line segments
% - Choose a mid-point in the data and compute the two best fit lines to
%   the data, with the constraint that the lines are perpendicular
% - Modify the inclusion of points in the two sets until the best fit is
%   found
% Modified Algorithm:
% - Determine the farthest two data points from each other
% - Use the adjacent data points in the set to determine the best fit lines
% - Work toward the intersection, point-by-point, adding points to
%   whichever set/line yields a lower residual
% NOTES:
% - Neither algorithm handles points on the "third" or "fourth" side of the
%   rectangle
% - Both algorithms seem computationally expensive, though the first
%   algorithm is closer to a bisection search and may be better.


% Determine the best fit (need to add rectangles and criteria for
% irregular)
for i_patch = 1:NumPatches
    patchArray(idxVec(i_patch)).primitive = 'circular';
end

% Plot, if requested by providing a figure handle
if flag_plot_prim
    figure(fig_num);
    
    % Check to see if hold is already on. If it is not, set a flag to turn it
    % off after this function is over so it doesn't affect future plotting
    flag_shut_hold_off = 0;
    if ~ishold
        flag_shut_hold_off = 1;
        hold on
    end
    
    % Plot the optimized circles
    theta = 0:2:360;
    hprim = zeros(NumPatches,1);
    for i_patch= 1:NumPatches
        hprim(i_patch) = plot(patchArray(idxVec(i_patch)).primparams(3)*cosd(theta)+patchArray(idxVec(i_patch)).primparams(1),...
            patchArray(idxVec(i_patch)).primparams(3)*sind(theta)+patchArray(idxVec(i_patch)).primparams(2),'k-.');
    end
    
    varargout{1} = hprim;
    
    % Shut the hold off?
    if flag_shut_hold_off
        hold off;
    end
end
