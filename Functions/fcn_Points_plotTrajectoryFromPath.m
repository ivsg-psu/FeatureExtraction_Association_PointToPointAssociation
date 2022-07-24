function fcn_Points_plotTrajectoryFromPath(path, fig_num)
% fcn_Points_plotLaneMarkers
% Adds elevation to the 2-dimensional 'path' based on the nearest neighbors
% in a 3-dimesional 'reference_elevated_path'
% 
% FORMAT:
%
%      fcn_Points_plotLaneMarkers(path,(fig_num))
%
% INPUTS:
%
%      path: a Nx2 vector of [X Y] path points, where N is the number of
%      points on the path, N >= 2.
%
%      (optional inputs)
%
%      figure_number: figure number where results are plotted
%
% OUTPUTS:
%
%      elevated_path: a Nx3 vector of [X Y Z] path points, where N is the 
%      number of points on the path, N >= 2.
%
% EXAMPLES:
% 
% See the script: script_test_fcn_Points_plotLaneMarkers for a full test 
% suite.
%
% DEPENDENCIES:
%
%
% This function was written on 2022_07_19 by Shashank Garikipati
% Questions or comments? sbrennan@psu.edu

% Revision history:
%      2022_07_19


flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _       
%  |_   _|                 | |      
%    | |  _ __  _ __  _   _| |_ ___ 
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |                  
%              |_| 
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Are the input vectors the right shape?
if 1 == flag_check_inputs
    % Are there the right number of inputs?
    if 2 ~= nargin
        error('Incorrect number of input arguments')
    end
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _       
%  |  \/  |     (_)      
%  | \  / | __ _ _ _ __  
%  | |\/| |/ _` | | '_ \ 
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LaneMarkers = 2;
pathToTraversal = fcn_Path_convertPathToTraversalStructure(path)

figure(fig_num);

fcn_Path_fillOffsetTraversalsAboutTraversal(pathToTraversal,[1;-1],fig_num);
%% Plot the results (for debugging)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _                 
%  |  __ \     | |                
%  | |  | | ___| |__  _   _  __ _ 
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % Ends the function