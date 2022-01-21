function h = fcn_Dataset_plotSetsXY(datasets,varargin)
% fcn_Dataset_plotSetsXY
% Plots the XY positions of all datasets existing in a data structure
%
% FORMAT: 
%
%       h = fcn_Dataset_plotSetsXY(datasets,{fig_num})
%
% INPUTS:
%
%      data: a structure containing subfields of X and Y coordinates in the
%      following form
%           data.set{i_set}.X
%           data.set{i_set}.Y
%      Note that i_path denotes an array of data sets. Each set will be
%      plotted separately.
%
% OUTPUTS:
%
%      h: a handle to the resulting figure
%
% DEPENDENCIES:
%
%      fcn_DataSet_checkInputsToFunctions
%
% EXAMPLES:
%      
%       See the script: script_test_fcn_Dataset_plotSetsXY.m for a full test
%       suite. 
%
% This function was adapted on 2022_01_12 by C. Beal from S. Brennan's
% fcn_Path_plotTraversalsXY
% Questions or comments? sbrennan@psu.edu 

% Revision history:
%     2022_01_12 
%     -- adapted for XY datasets without path ordering

flag_do_debug = 0; % Flag to plot the results for debugging
flag_this_is_a_new_figure = 1; % Flag to check to see if this is a new figure
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
    if nargin < 1 || nargin > 2
        error('Incorrect number of input arguments')
    end
    
    % Check the data input
    % fcn_Path_checkInputsToFunctions(traversals, 'traversals');

end

% Does user want to show the plots?
if 2 == nargin
    fig_num = varargin{1};
    figure(fig_num);
    flag_this_is_a_new_figure = 0;
else    
    fig = figure;
    fig_num = fig.Number;
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

NumSets = length(datasets.traversal);
h = zeros(NumTraversals,1);
for i_path= 1:NumTraversals
    h(i_path) = plot(datasets.traversal{i_path}.X,datasets.traversal{i_path}.Y,'-o');
end

% Shut the hold off?
if flag_shut_hold_off
    hold off;
end

% Add labels? 
if flag_this_is_a_new_figure == 1
    title('X vs Y')
    xlabel('X [m]')
    ylabel('Y [m]')
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

