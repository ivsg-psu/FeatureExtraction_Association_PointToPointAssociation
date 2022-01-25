function [h, hpts] = fcn_Dataset_plotPatch(patchArray,varargin)
% fcn_Dataset_plotPatch
% Plots a visual representation of the objects in a patch structure array
%
% FORMAT: 
%
%       [h,hpts] = fcn_Dataset_plotPatch(patchArray,{fig_num},{indices})
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
%      (OPTIONAL INPUTS)
%      fig_num: a figure number to plot into
%      indices: a vector of indices indicating which patch objects to plot
%
% OUTPUTS:
%
%      h: a vector of handles to the resulting patch objects plotted in
%         order of the indices (if provided)
%      hpts: a vector of handles to the resulting points plotted on the
%            patch in order of the indices (if provided)
%
% DEPENDENCIES:
%
%      ## NOT CURRENTLY USED: fcn_DataSet_checkInputsToFunctions
%
% EXAMPLES:
%      
%       See the script: script_test_fcn_Dataset_plotPatch.m for a full test
%       suite. 
%
% This function was written by C. Beal
% Questions or comments? cbeal@bucknell.edu 

% Revision history:
%     2022_01_25
%     -- wrote the code

flag_do_debug = 0; % Flag to plot the results for debugging
flag_this_is_a_new_figure = 1; % Flag to check to see if this is a new figure
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
    flag_this_is_a_new_figure = 0;
else    
    fig = figure;
    fig_num = fig.Number;
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
 
figure(fig_num);

% Check to see if hold is already on. If it is not, set a flag to turn it
% off after this function is over so it doesn't affect future plotting
flag_shut_hold_off = 0;
if ~ishold
    flag_shut_hold_off = 1;
    hold on
end

NumPatches = length(idxVec);
h = zeros(NumPatches,1);
hpts = zeros(NumPatches,1);
for i_patch= 1:NumPatches
    h(i_patch) = patch(patchArray(idxVec(i_patch)).pointsX,patchArray(idxVec(i_patch)).pointsY,patchArray(idxVec(i_patch)).color);
    hpts(i_patch) = plot(patchArray(idxVec(i_patch)).pointsX,patchArray(idxVec(i_patch)).pointsY,'k*');
end

% Shut the hold off?
if flag_shut_hold_off
    hold off;
end

% Add labels? 
if flag_this_is_a_new_figure == 1
    axis equal
    grid on
    title('X vs Y')
    xlabel('X [m]')
    ylabel('Y [m]')
end

