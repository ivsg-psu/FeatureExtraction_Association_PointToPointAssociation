function [patchArray, varargout] = fcn_Dataset_determineAABB(patchArray,varargin)
% fcn_Dataset_determineAABB
% Determine the axis-aligned bounding box for a patch structure or array of
% patch structures
%
% FORMAT:
%
%       patchArray = fcn_Dataset_determineAABB(patchArray,{fig_num},{indices})
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
%      hbb: a vector of handles to the bounding boxes in order of the 
%           indices (if provided)
%
% DEPENDENCIES:
%
%      ## NOT CURRENTLY USED: fcn_DataSet_checkInputsToFunctions
%
% EXAMPLES:
%
%       See the script: script_test_fcn_Dataset_determineAABB.m for a full test
%       suite.
%
% This function was written by C. Beal
% Questions or comments? cbeal@bucknell.edu

% Revision history:
%     2022_01_25
%     -- wrote the code

flag_do_debug = 0; % Flag to add debugging code
flag_check_inputs = 1; % Flag to perform input checking
flag_plot_aabb = 0; % Flag to plot the axis-aligned bounding boxes

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
    flag_plot_aabb = 1;
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
for i_patch = 1:NumPatches
    patchArray(idxVec(i_patch)).aabb = [min(patchArray(idxVec(i_patch)).pointsX) min(patchArray(idxVec(i_patch)).pointsY)...
        max(patchArray(idxVec(i_patch)).pointsX) max(patchArray(idxVec(i_patch)).pointsY)];
end

% Plot, if requested by providing a figure handle
if flag_plot_aabb
    figure(fig_num);
    
    % Check to see if hold is already on. If it is not, set a flag to turn it
    % off after this function is over so it doesn't affect future plotting
    flag_shut_hold_off = 0;
    if ~ishold
        flag_shut_hold_off = 1;
        hold on
    end
    
    hbb = zeros(NumPatches,1);
    for i_patch= 1:NumPatches
        hbb(i_patch) = plot(patchArray(idxVec(i_patch)).aabb([1 3 3 1 1]),patchArray(idxVec(i_patch)).aabb([2 2 4 4 2]),'k-','linewidth',2);
    end
    
    varargout{1} = hbb;
    
    % Shut the hold off?
    if flag_shut_hold_off
        hold off;
    end
end
