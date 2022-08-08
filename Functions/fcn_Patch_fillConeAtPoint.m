function [conePatch] = fcn_Patch_fillConeAtPoint(xylocation,orientation)
% fcn_Patch_fillBarrelAtPoint
% Plots a 2D barrel at the user's desired location with the user's desired
% orientation
%
% FORMAT:
%
%       samplePatches = fcn_Patch_fillConeAtPoint(location,orientation)
%
% INPUTS:
%      numberOfObstables: 1x1 scalar representin number of obstacles
%      xylocation: [1x2] matrix of [x y] points, representing the location of
%      the cone
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
%      2022_08_08 - George Bodenschatz
%      -- Added visual detail to cones

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
conePatch = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});

white = [1 1 1];
orange = [1 .65 0]
black = [0 0 0]

for i = 1:size(xylocation,1)
    conePatch(i).color = orange;
    conePatch(i+size(xylocation,1)).color = white;
    conePatch(i+2*size(xylocation,1)).color = orange;
    conePatch(i+3*size(xylocation,1)).color = white;
    conePatch(i+4*size(xylocation,1)).color = orange;
    conePatch(i+5*size(xylocation,1)).color = black;
    r = 0.18;
    s = r*2
    
      p1_x = -s/2;
    p1_y = s/2;

    p2_x = -s/2;
    p2_y = -s/2

    p3_x = s/2;
    p3_y = -s/2

    p4_x = s/2;
    p4_y = s/2;
    
    t = linspace(0,2*pi);
    xycirclepoints = [cos(t);sin(t)];
    
    rot_mat = [cos(orientation(i)) -sin(orientation(i)); sin(orientation(i)) cos(orientation(i))];
    
    xypoints = rot_mat*[p1_x p2_x p3_x p4_x; p1_y p2_y p3_y p4_y];
    
    conePatch(i).pointsX = [xypoints(1,:)]' + xylocation(i,1);
    conePatch(i).pointsY = [xypoints(2,:)]' + xylocation(i,2);


    conePatch(i+size(xylocation,1)).pointsX = [(r-.03)*xycirclepoints(1,:)]' + xylocation(i,1);
    conePatch(i+size(xylocation,1)).pointsY = [(r-.03)*xycirclepoints(2,:)]' + xylocation(i,2);
    conePatch(i+2*size(xylocation,1)).pointsX = [(r-.06)*xycirclepoints(1,:)]' + xylocation(i,1);
    conePatch(i+2*size(xylocation,1)).pointsY = [(r-.06)*xycirclepoints(2,:)]' + xylocation(i,2);
    conePatch(i+3*size(xylocation,1)).pointsX = [(r-.09)*xycirclepoints(1,:)]' + xylocation(i,1);
    conePatch(i+3*size(xylocation,1)).pointsY = [(r-.09)*xycirclepoints(2,:)]' + xylocation(i,2);
    conePatch(i+4*size(xylocation,1)).pointsX = [(r-.12)*xycirclepoints(1,:)]' + xylocation(i,1);
    conePatch(i+4*size(xylocation,1)).pointsY = [(r-.12)*xycirclepoints(2,:)]' + xylocation(i,2);
    conePatch(i+5*size(xylocation,1)).pointsX = [(r-.15)*xycirclepoints(1,:)]' + xylocation(i,1);
    conePatch(i+5*size(xylocation,1)).pointsY = [(r-.15)*xycirclepoints(2,:)]' + xylocation(i,2);

    
end

