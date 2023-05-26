function datasets_out = fcn_Points_addRadialNoise(datasets,R)
% fcn_Points_addRadialNoise(datasets,R)
% Adds radial noise to the datasets using R (radius of max noise)
%
% FORMAT: 
%
%       h = fcn_Points_addRadialNoise(datasets,R)
%
% INPUTS:
%
%      datasets: a structure array containing subfields of X and Y 
%      coordinates in the following form:
%           datasets{i_set}.X
%           datasets{i_set}.Y
%      Note that i_set addresses a particular data set structure. Each set
%      will be modified separately.
%      R: radius of max noise
%
%                   OPTIONAL INPUTS:
%
%      fig_num: figure number
%
% OUTPUTS:
%
%      datasets_out: the input datasets will be modified and returned
%
% DEPENDENCIES:
%
%      fcn_Points_checkInputsToFunctions
%
% EXAMPLES:
%      
%       See the script: script_test_fcn_Points_adjustPointSetStatistics.m
%       for a full test suite. 
%
% This function was written on 2022_07_18 by Shashank Garikipati
% Questions or comments?

% Revision history:
%     2022_07_18
%     -- wrote the code

flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end


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

if flag_check_inputs == 1
    % Are there the right number of inputs?
    if 2 < nargin || nargin > 3
        error('Incorrect number of input arguments')
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
 
% Iterate through each of the datasets and add the radial noise

x = zeros([size(datasets,2) 1]);
y = zeros([size(datasets,2) 1]);
datasets_out = {zeros([size(datasets,2) 2])};
for i_set= 1:size(datasets,2)
    for i = 1:size(datasets{i_set},1)
        rand_radius = R*rand;
        rand_angle = 2*pi*rand;
        x_data = datasets{i_set}(i,1);
        y_data = datasets{i_set}(i,2);
        x(i,1) = rand_radius*cos(rand_angle) + x_data;
        y(i,1) = rand_radius*sin(rand_angle) + y_data;
    end
    datasets_out{i_set} = [x y];
end


%% Any debugging?
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
if flag_do_debug
    % Nothing in here yet
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
end

