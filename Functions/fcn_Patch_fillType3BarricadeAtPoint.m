function [type3Barricade] = fcn_Patch_fillConeAtPoint(xylocation,orientation)
% fcn_Patch_fillBarrelAtPoint
% Plots a 2D barrel at the user's desired location with the user's desired
% orientation
%
% FORMAT:
%
%       samplePatches = fcn_Patch_fillBarrelAtPoint(location,orientation)
%
% INPUTS:
%      xylocation: [1x2] matrix of [x y] points, representing the location of
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
%       script_test_fcn_Patch_fillType3BarricadeAtPoint.m for a full
%       test suite.
%
% This function was written on 2022_07_15 by Shashank Garikipati

% Revision history:
%      2020_07_15
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
type3Barricade = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});

for i = 1:size(xylocation,1)
    type3Barricade(i).color = [.5 .5 .5];
    type3Barricade(i+size(xylocation,1)).color = [.5 .5 .5]
    type3Barricade(i+2*size(xylocation,1)).color = [1 .65 0]
    
    height = 0.1;
    width = 1.2;
    standheight = 1;
    standwidth = .05;
    
    p1_x = -width/2;
    p1_y = -height/2;

    p2_x = -width/2;
    p2_y = height/2 ;

    p3_x = width/2;
    p3_y = height/2;

    p4_x = width/2;
    p4_y = -height/2;
    
    stand1p1_x = width/2;
    stand1p1_y = standheight/2;
    
    stand1p2_x = width/2;
    stand1p2_y = -standheight/2;
    
    stand1p3_x = width/2 - standwidth;
    stand1p3_y = -standheight/2;
    
    stand1p4_x = width/2 - standwidth;
    stand1p4_y = standheight/2;
    
    stand2p1_x = -width/2;
    stand2p1_y = standheight/2;
    
    stand2p2_x = -width/2;
    stand2p2_y = -standheight/2;
    
    stand2p3_x = -(width/2 - standwidth);
    stand2p3_y = -standheight/2;
    
    stand2p4_x = -(width/2 - standwidth);
    stand2p4_y = standheight/2;
    
    rot_mat = [cos(orientation(i)) -sin(orientation(i)); sin(orientation(i)) cos(orientation(i))];
    
    xypoints = rot_mat*[p1_x p2_x p3_x p4_x; p1_y p2_y p3_y p4_y];
    xystand1points = rot_mat*[stand1p1_x stand1p2_x stand1p3_x stand1p4_x; stand1p1_y stand1p2_y stand1p3_y stand1p4_y];
    xystand2points = rot_mat*[stand2p1_x stand2p2_x stand2p3_x stand2p4_x; stand2p1_y stand2p2_y stand2p3_y stand2p4_y];

    type3Barricade(i).pointsX = [xystand2points(1,:)]' + xylocation(i,1);
    type3Barricade(i).pointsY = [xystand2points(2,:)]' + xylocation(i,2);
    type3Barricade(i+size(xylocation,1)).pointsX = [xystand1points(1,:)]' + xylocation(i,1);
    type3Barricade(i+size(xylocation,1)).pointsY = [xystand1points(2,:)]' + xylocation(i,2);
    type3Barricade(i+2*size(xylocation,1)).pointsX = [xypoints(1,:)]' + xylocation(i,1);
    type3Barricade(i+2*size(xylocation,1)).pointsY = [xypoints(2,:)]' + xylocation(i,2);
end

