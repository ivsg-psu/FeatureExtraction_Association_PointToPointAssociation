function h = fcn_Points_plotSetsXY(datasets,varargin)
% fcn_Points_plotSetsXY
% Plots the XY positions of all datasets existing in a data structure
%
% FORMAT: 
%
%       h = fcn_Points_plotSetsXY(datasets,{fig_num})
%
% INPUTS:
%
%      datasets: a structure array containing subfields of X and Y 
%      coordinates in the following form:
%           datasets{i_set}.X
%           datasets{i_set}.Y
%      Note that i_set denotes a data set structure. Each set will be
%      plotted separately.
%
% OUTPUTS:
%
%      h: a handle to the resulting figure
%
% DEPENDENCIES:
%
%      fcn_Points_checkInputsToFunctions
%
% EXAMPLES:
%      
%       See the script: script_test_fcn_Points_plotSetsXY.m for a full test
%       suite. 
%
% This function was adapted on 2022_01_12 by C. Beal from S. Brennan's
% fcn_Path_plotTraversalsXY
% Questions or comments? sbrennan@psu.edu 

% Revision history:
% 2022_01_12 
% -- adapted for XY datasets without path ordering
% 2023_05_29
% -- fixed latent bug if user enters empty figure number

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
    narginchk(1,2);
    
    % Check the data input
    % fcn_Path_checkInputsToFunctions(traversals, 'traversals');

end

% Does user want to show the plots?
flag_this_is_a_new_figure = 1; % Flag to check to see if this is a new figure
if 2 == nargin
    temp = varargin{1};
    if ~isempty(temp)
        fig_num = temp;
        figure(fig_num);
        flag_this_is_a_new_figure = 0;
    end
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

NumSets = length(datasets);
h = zeros(NumSets,1);
for i_set= 1:NumSets
    try
    h(i_set) = plot(datasets{i_set}(:,1),datasets{i_set}(:,2),'o');
    catch
        disp('stop here');
    end
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

