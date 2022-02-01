function patchStruct = fcn_Patch_fillPatchArrayViaUserInputs(varargin)
% fcn_Patch_fillPatchArrayViaUserInputs
% A function for the user to click on the figure to generate XY data.
% Points are collected and plotted until the user double clicks. If the
% user right-clicks anywhere in the plot, the last point is deleted. Once
% the user double-clicks, the results are output from the function.
%
% FORMAT:
%
%      patchStruct = fcn_Patch_fillPatchArrayViaUserInputs({fig_num})
%
% INPUTS:
%
%      (OPTIONAL INPUTS)
%      fig_num: an integer specifying which figure to use
%
% OUTPUTS:
%
%      patchStruct: structure array patch objects that the user generated
%                  via clicking in a map and selecting properties
%
% EXAMPLES:
%
%      % BASIC examples
%      patchStruct = fcn_Patch_fillPatchArrayViaUserInputs
%      patchStruct = fcn_Patch_fillPatchArrayViaUserInputs(1)
%
% See the script: script_test_fcn_Patch_fillPatchArrayViaUserInputs
% for a full test suite.
%
% This function was adapted on 2022_01_28 by C. Beal from S. Brennan's
% fcn_Path_fillPathViaUserInputs function from the PSU IVSG path library.
% Questions or comments? sbrennan@psu.edu

% Revision history:
%      2022_01_28
%      -- adapted for patch objects in a structure array

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
    fig_num = varargin{1};
    callback_details = varargin{2};
    callback_type = callback_details.EventName;
elseif 0 == nargin
    fig_num = figure;
else
    fig_num = varargin{1};
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
    % Make a new patch structure
    patchStruct = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});
    % Fill out the patch id with a new identifier
    [~,patchStruct(1).id] = fileparts(tempname);
    % Define the patch as having an irregular primitive
    patchStruct(1).primitive = 'irregular';
    % Set the patch color
    patchStruct(1).color = uisetcolor();
    % Make a new figure, initializing all data and handles within
    fcn_Patch_fillPatchArrayViaUserInputs_startPlot(patchStruct(1),fig_num);
    UserData = get(gcf,'UserData');
    % This loop will continue to run until the callback for a double-click
    % is handled to end point addition
    while UserData.flag_is_done == 0
        % Wait for the figure to be done
        pause(0.15);
        UserData = get(gcf,'UserData');
    end
    % Set the patch point data 
    patchStruct(1) = UserData.patch;
    % When this portion of the code finishes, the function is complete
    % and the patch structure will be returned
    
else
    switch callback_type
        case {'WindowMouseMotion'}
            UserData = get(gcf,'UserData');
            prompt = UserData.title_header;
            
            C = get (gca, 'CurrentPoint');
            title_string = sprintf('%s, (X,Y) = %.1f, %.1f',prompt,C(1,1),C(1,2));
            title(gca, title_string);
        case {'WindowMousePress'}
            % Mouse was clicked - find the location
            mousePos = get (gca, 'CurrentPoint');
            
            % Find the type of click
            source = callback_details.Source;
            mouseType = get(source, 'SelectionType');
            
            switch mouseType
                case {'normal'}  % Typical left-click - adds a point (even if a double-click is also processed)
                    if flag_do_debug
                        fprintf(1,'Left clicked -  X: %.1f, y: %.1f\n',mousePos(1),mousePos(2));
                    end
                    
                    % Grab the data out of the figure so we can update the data and the plot
                    UserData = get(gcf,'UserData');
                    
                    % Insert the new point into the patch structure
                    UserData.patch = fcn_Patch_insertPoints(UserData.patch,mousePos(1,1:2));
                    
                    % Update the plot from the patch structure
                    delete(UserData.h_patch);
                    UserData.h_patch = patch(UserData.patch.pointsX,UserData.patch.pointsY,UserData.patch.color);
                    
                    % Save the results so they can be accessed by the
                    % various callbacks and force a plot update
                    set(gcf,'UserData',UserData);
                    drawnow
                                        
                case {'alt'}  % Typical right-click - subtracts a point
                    if flag_do_debug
                        fprintf(1,'Right clicked -  X: %.1f, y: %.1f\n',mousePos(1),mousePos(2));
                    end
                    
                    % Grab the data out of the figure so we can update the plot
                    UserData = get(gcf,'UserData');
                    
                    % Mouse was clicked - find the location
                    mousePos = get(gca,'CurrentPoint');
            
                    % Find the closest point to the mouse click
                    idx = knnsearch([UserData.patch.pointsX UserData.patch.pointsY],mousePos(1,1:2));
                    
                    % Remove the closest point and shift the data back
                    % Handle the case where the index is at the end
                    if idx == length(UserData.patch.pointsX)
                        UserData.patch.pointsX = UserData.patch.pointsX(1:end-1);
                        UserData.patch.pointsY = UserData.patch.pointsY(1:end-1);
                    % Handle the case where the index is at the start
                    elseif 1 == idx
                        UserData.patch.pointsX = UserData.patch.pointsX(2:end);
                        UserData.patch.pointsY = UserData.patch.pointsY(2:end);
                    % Handle the case where the index is in the midst of
                    % the vector
                    else
                        UserData.patch.pointsX = [UserData.patch.pointsX(1:idx-1); UserData.patch.pointsX(idx+1:end)];
                        UserData.patch.pointsY = [UserData.patch.pointsY(1:idx-1); UserData.patch.pointsY(idx+1:end)];
                    end
                    
                    % Update the plot from the patch structure
                    delete(UserData.h_patch);
                    UserData.h_patch = patch(UserData.patch.pointsX,UserData.patch.pointsY,UserData.patch.color);
                    
                    % Save the results so they can be accessed by the
                    % various callbacks and force a plot update
                    set(gcf,'UserData',UserData);
                    drawnow
                    
                case {'open'} % Typical double click - ends the routine
                    if flag_do_debug
                        fprintf(1,'Double clicked -  X: %.1f, y: %.1f\n',mousePos(1),mousePos(2));
                    end
                    
                    % Shut off callbacks
                    set(gcf,'WindowButtonMotionFcn',{}, ...
                        'WindowButtonDownFcn',{});
                    
                    % Grab the data out of the figure so we can update it
                    UserData = get(gcf,'UserData');
                    
                    % Set the flag to update the preceding callbacks that
                    % the process is complete
                    UserData.flag_is_done = 1;
                    
                    % Save the results so they can be accessed by the
                    % various callbacks and force a plot update
                    set(gcf,'UserData',UserData);
                    
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

function fcn_Patch_fillPatchArrayViaUserInputs_startPlot(patchStruct,fig_num)
% Decide the number of points to use (maximum), initialize data and values
num_points = 1000;
data = nan*ones(num_points,2);
next_point = 1;

% Set up the figure
current_fig = figure(fig_num);

% Plot empty patch
h_patch = patch([],[],patchStruct.color);
xlim([0 100]);
ylim([0 100]);

% Grab the data out of the figure so we can update it
UserData = get(gcf,'UserData');

% Fill in a UserData structure
UserData.num_points = num_points;
UserData.data = data;
UserData.next_point = next_point;
UserData.h_patch = h_patch;
UserData.flag_is_done = 0;
UserData.patch = patchStruct;

% Check to see that a title has been provided
if ~isfield(UserData,'title_header')
    UserData.title_header = sprintf('Enter XY point data by clicking');
end

% Save Userdata to the figure
set(current_fig,'UserData',UserData);

% Make sure the figure updates with the limits, title, axes labels, etc.
drawnow

% Save user data into the current figure, associating movement and clicking
% functions to specific functions
set(current_fig, 'WindowButtonMotionFcn',...
    @fcn_Patch_fillPatchArrayViaUserInputs, ...
    'WindowButtonDownFcn',@fcn_Patch_fillPatchArrayViaUserInputs);
end