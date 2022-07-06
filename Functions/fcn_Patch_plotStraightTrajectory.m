function fcn_Patch_plotStraightTrajectory(t_h,x0,vehicle,varargin)
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
%      x0: [4x1] vector of [p0',h0,vx]
%         vx: [1x1] scalar, representing longitudinal speed
%         p0: [1x1] scalar, representing initial position of the vehicle
%         h0: [1x1] scalar, representing initial heading of the vehicle
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
    if nargin < 3 || nargin > 4
        error('Incorrect number of input arguments')
    end
    
end

% Did the user provide a figure number?
if 3 < nargin
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

p0 = x0(1:2)';
h0 = x0(3);
vx = x0(4);

plot(p0(1),p0(2),'k*')

del_t = 0.01;
t = (0:del_t:t_h)';
d = vx*t;
N = length(t);
pv = zeros(N,2);
pv(:,1) = d*cos(h0) + p0(1);
pv(:,2) = d*sin(h0) + p0(2);
% And plot the trajectory
plot(pv(:,1),pv(:,2),'k-.')

% plot outer limits of vehicle trajectory
plot(pv(:,1)-vehicle.dr*cos(h0)-vehicle.w/2*sin(h0),pv(:,2)-vehicle.dr*sin(h0)+vehicle.w/2*cos(h0),'b');
plot(pv(:,1)-vehicle.dr*cos(h0)+vehicle.w/2*sin(h0),pv(:,2)-vehicle.dr*sin(h0)-vehicle.w/2*cos(h0),'r');



