function [barrelPatch] = fcn_Patch_fillBarrelAtPoint(xylocation,orientation)
% fcn_Patch_fillBarrelAtPoint
% Plots a 2D barrel at the user's desired location with the user's desired
% orientation
%
% FORMAT:
%
%       samplePatches = fcn_Patch_fillBarrelAtPoint(location,orientation)
%
% INPUTS:
%      numberOfObstables: 1x1 scalar representin number of obstacles
%      xylocation: [Nx2] matrix of [x y] points, representing the location of
%      the barrel
%      orientation: -----------,representing the orientation of the barrel
%
% OUTPUTS:
%
%      barrelPatch: an cell array of datasets
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_Patch_fillBarrelAtPoint.m for a full
%       test suite.
%
% This function was written on 2022_07_05 by Shashank Garikipati

% Revision history:
%      2020_01_28
%      -- wrote the code

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
    if  nargin > 2 || nargin < 2
        error('Incorrect number of input arguments')
    end
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

%% create new cell structure for barrel information
barrelPatch = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});

for i = 1:size(xylocation,1)
    barrelPatch(i).color = [1 .65 0];
    r = 0.2;
    
    p1_x = (-r/(sqrt(2)+1));
    p1_y = -r;
    
    p2_x = -r;
    p2_y = (-r/(sqrt(2)+1));
    
    p3_x = -r;
    p3_y = (r/(sqrt(2)+1));
    
    p4_x = (-r/(sqrt(2)+1));
    p4_y = r;
    
    p5_x = (r/(sqrt(2)+1));
    p5_y = r;
    
    p6_x = r;
    p6_y = (r/(sqrt(2)+1));
    
    p7_x = r;
    p7_y = -(r/(sqrt(2)+1));
    
    p8_x = (r/(sqrt(2)+1));
    p8_y = -r;
    
    rot_mat = [cos(orientation(i)) -sin(orientation(i)); sin(orientation(i)) cos(orientation(i))];
    
    xypoints = rot_mat*[p1_x p2_x p3_x p4_x p5_x p6_x p7_x p8_x; p1_y p2_y p3_y p4_y p5_y p6_y p7_y p8_y];
    
    barrelPatch(i).pointsX = [xypoints(1,:)]' + xylocation(i,1);
    barrelPatch(i).pointsY = [xypoints(2,:)]' + xylocation(i,2);
end

