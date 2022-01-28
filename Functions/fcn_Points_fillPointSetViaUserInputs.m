function dataXY = fcn_Points_fillPointSetViaUserInputs(fig_num,varargin)
% fcn_Points_fillPointSetViaUserInputs
% A function for the user to click on the figure to generate XY data.
% Points are collected and plotted until the user double clicks. If the
% user right-clicks anywhere in the plot, the last point is deleted. Once
% the user double-clicks, the results are output from the function.
%
% FORMAT:
%
%      dataXY = fcn_Points_fillPointSetViaUserInputs({fig_num})
%
% INPUTS:
%
%      (OPTIONAL INPUTS)
%      fig_num: an integer specifying which figure to use
%
% OUTPUTS:
%
%      dataXY: matrix (Nx2) representing the X and Y points that the user
%      clicked on the map
%
% EXAMPLES:
%
%      % BASIC examples
%      dataXY = fcn_Points_fillPointSetViaUserInputs
%      dataXY = fcn_Points_fillPointSetViaUserInputs(1)
%
% See the script: script_test_fcn_Points_fillPointSetViaUserInputs
% for a full test suite.
%
% This function was adapted on 2022_01_12 by C. Beal from S. Brennan's
% fcn_Path_fillPathViaUserInputs function from the PSU IVSG path library.
% Questions or comments? sbrennan@psu.edu

% Revision history:
%      2022_01_12
%      -- adapted for XY datasets without path ordering

flag_do_debug = 0; % Flag to plot the results for debugging
flag_make_figure = 0;
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    flag_make_figure = 1;
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
    if nargin > 2
        error('Incorrect number of input arguments')
    end 
end

callback_type = '';
% Is this a call-back?
if 2 == nargin
    callback_details = varargin{1};
    callback_type = callback_details.EventName;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(callback_type)
    % Make a new figure, initializing all data and handles within
    fcn_Points_fillPointSetViaUserInputs_startPlot(fig_num);
    UserData = get(gcf,'UserData');
    while UserData.flag_is_done == 0
        % Wait for the figure to be done
        pause(0.15);
        UserData = get(gcf,'UserData');
    end
    current_point = UserData.next_point;
    dataXY = UserData.data(1:current_point-1,:);
    
else
    switch callback_type
        case {'WindowMouseMotion'}
            UserData = get(gcf,'UserData');
            prompt = UserData.title_header;
            
            % Mouse was moved
            C = get (gca, 'CurrentPoint');
            title_string = sprintf('%s, (X,Y) = %.1f, %.1f',prompt,C(1,1),C(1,2));
            title(gca, title_string);
        case {'WindowMousePress'}
            % Mouse was clicked - find the location
            mousePos = get (gca, 'CurrentPoint');
            
            % Find the type of click
            source   = callback_details.Source;
            mouseType = get(source, 'SelectionType');
            
            switch mouseType
                case {'normal'}  % Typical left-click - adds a point
                    if flag_do_debug
                        fprintf(1,'Left clicked -  X: %.1f, y: %.1f\n',mousePos(1),mousePos(2));
                    end
                    
                    % Grab the data out of the figure so we can update the plot
                    UserData = get(gcf,'UserData');
                    num_points = UserData.num_points;
                    data = UserData.data;
                    
                    % Fill in the current point
                    current_point = UserData.next_point;
                    data(current_point,1) = mousePos(1,1);
                    data(current_point,2) = mousePos(1,2);
                    
                    % Find the point after the current point, and shift data down if "full"
                    next_point = current_point+1;
                    if next_point == (num_points+1)
                        data(1:end-1,:) = data(2:end,:);
                        next_point = num_points;
                    end
                    
                    % Update the plot
                    set(UserData.h_plot,'XData',data(:,1),'YData',data(:,2));
                    
                    %UserData.num_points = num_points;
                    UserData.data = data;
                    UserData.next_point = next_point;
                    
                    % Save the results
                    set(gcf,'UserData',UserData);
                case {'alt'}  % Typical right-click - subtracts a point
                    if flag_do_debug
                        fprintf(1,'Right clicked -  X: %.1f, y: %.1f\n',mousePos(1),mousePos(2));
                    end
                    
                    % Grab the data out of the figure so we can update the plot
                    UserData = get(gcf,'UserData');
                    num_points = UserData.num_points;
                    data = UserData.data;
                    
                    % Go backwards
                    UserData.next_point = max(UserData.next_point - 1,1);
                    
                    % Fill in the current point
                    current_point = UserData.next_point;
                    data(current_point,1) = nan;
                    data(current_point,2) = nan;
                    
                    % Update the plot
                    set(UserData.h_plot,'XData',data(:,1),'YData',data(:,2));
                    
                    % Update the data
                    UserData.data = data;
                    
                    % Save the results
                    set(gcf,'UserData',UserData);
                    
                    
                case {'open'} % Typical double click - ends the routine
                    if flag_do_debug
                        fprintf(1,'Double clicked -  X: %.1f, y: %.1f\n',mousePos(1),mousePos(2));
                    end
                    
                    % Shut off callbacks
                    set(gcf, ...
                        'WindowButtonMotionFcn',{}, ...
                        'WindowButtonDownFcn',{});
                    
                    % Grab the data out of the figure so we can update it
                    UserData = get(gcf,'UserData');
                    
                    % Set the flag that we are done
                    UserData.flag_is_done = 1;
                    
                    % Save the results
                    set(gcf,'UserData',UserData);
                    
                    if flag_do_debug
                        disp('DONE');
                    end
                    
                otherwise
                    fprintf(1,'Unknown mouse click type: %s\n',mouseType);
            end
            
            
            
        otherwise
            fprintf(1,'unknown callback type: %s\n',callback_type);
    end
end

%%
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
if flag_make_figure
    %     figure(fig_num);
    %     hold on;
    %     grid on;
    %     plot(dataXY(:,1),dataXY(:,2),'r.','Markersize',20);
    %     plot(dataXY(:,1),dataXY(:,2),'r-','Linewidth',3);
    
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end
end

function fcn_Points_fillPointSetViaUserInputs_startPlot(fig_num)
% Decide the number of points to use (maximum), initialize data and values
num_points = 1000;
data = nan*ones(num_points,2);
next_point = 1;
dataXY = [0 0; 0 1];

% Set up the figure
current_fig = figure(fig_num);

% Plot empty data
h_plot = plot(data(:,1),data(:,2),'.','Markersize',20);
xlim([0 100]);
ylim([0 100]);

% Grab the data out of the figure so we can update it
UserData = get(gcf,'UserData');

% Fill in a UserData structure
UserData.num_points = num_points;
UserData.data = data;
UserData.next_point = next_point;
UserData.h_plot = h_plot;
UserData.flag_is_done = 0;

% Save Userdata to the figure
set(current_fig,'UserData',UserData);

% Save user data into the current figure, associating movement and clicking
% functions to specific functions
set(current_fig, 'WindowButtonMotionFcn',...
    @fcn_Points_fillPointSetViaUserInputs, ...
    'WindowButtonDownFcn',@fcn_Points_fillPointSetViaUserInputs);
end





