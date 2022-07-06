function fcn_Patch_plotStraightTrajectory(collFlags,collLoc,clearance,bodyCollLoc,vehicle,h0,varargin)
% fcn_Patch_plotPatch
% Plots a visual representation of the objects in a patch structure array
%
% FORMAT: 
%
%       fcn_Patch_plotCenterGravity(centerOfGravity)
%
% INPUTS:
%
%      tf: [1x1] scalar, representing time horizon to check
%      vx: [1x1] scalar, representing longitudinal speed
%      CG: [1x1] scalar, representing initial position of the vehicle
%      h0: [1x1] scalar, representing initial heading of the vehicle
%      vehicle: [1x1] cell struct with vehicle properties
%
%      (OPTIONAL INPUTS)
%      fig_num: a figure number to plot into
%
% OUTPUTS:
%
%
% DEPENDENCIES:
%
%
% EXAMPLES:
%      
%       See the script: script_test_fcn_Patch_plotStraightTragectory.m for a full test
%       suite. 
%
% This function was written by Shashank Garikipati

% Revision history:
%     2022_07_06
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
    if nargin < 6 || nargin > 7
        error('Incorrect number of input arguments')
    end
    
end

% Did the user provide a figure number?
if 6 < nargin
    fig_num = varargin{1};
    figure(fig_num);
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
grid on
hold on
axis equal

for collInd = 1:length(collFlags)
    plotOffset(1) = bodyCollLoc(collInd,1);
    if isnan(clearance(collInd))
        plotOffset(2) = bodyCollLoc(collInd,2);
    else
        plotOffset(2) = bodyCollLoc(collInd,2) + sign(bodyCollLoc(collInd,2))*clearance(collInd);
    end
    bodyLoc(1,1) = collLoc(collInd,1) + (vehicle.df - plotOffset(1))*cos(h0) - (vehicle.w/2-plotOffset(2))*sin(h0);
    bodyLoc(1,2) = collLoc(collInd,2) + (vehicle.df - plotOffset(1))*sin(h0) + (vehicle.w/2-plotOffset(2))*cos(h0);
    bodyLoc(2,1) = collLoc(collInd,1) + (vehicle.df - plotOffset(1))*cos(h0) - (-vehicle.w/2-plotOffset(2))*sin(h0);
    bodyLoc(2,2) = collLoc(collInd,2) + (vehicle.df - plotOffset(1))*sin(h0) + (-vehicle.w/2-plotOffset(2))*cos(h0);
    bodyLoc(3,1) = collLoc(collInd,1) + (-vehicle.dr - plotOffset(1))*cos(h0) - (-vehicle.w/2-plotOffset(2))*sin(h0);
    bodyLoc(3,2) = collLoc(collInd,2) + (-vehicle.dr - plotOffset(1))*sin(h0) + (-vehicle.w/2-plotOffset(2))*cos(h0);
    bodyLoc(4,1) = collLoc(collInd,1) + (-vehicle.dr - plotOffset(1))*cos(h0) - (vehicle.w/2-plotOffset(2))*sin(h0);
    bodyLoc(4,2) = collLoc(collInd,2) + (-vehicle.dr - plotOffset(1))*sin(h0) + (vehicle.w/2-plotOffset(2))*cos(h0);
    plot(bodyLoc([1:end 1],1),bodyLoc([1:end 1],2),'b-','linewidth',1);
    
    % Plot the collision or closest clearance point
    plot(collLoc(collInd,1),collLoc(collInd,2),'r*')
end



