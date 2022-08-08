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
%      2022_08_08 - George Bodenschatz
%      -- Added visual detail to barrels

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
    barrelPatch(i+size(xylocation,1)).color = [1 1 1];
    barrelPatch(i+2*size(xylocation,1)).color = [1 .65 0];
    barrelPatch(i+3*size(xylocation,1)).color = [1 1 1];
    barrelPatch(i+4*size(xylocation,1)).color = [1 .65 0];
    barrelPatch(i+5*size(xylocation,1)).color = [0 0 0];
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
    
        rectp1_x = 0.03;
    rectp1_y = 0.01;
    
    rectp2_x = -0.03;
    rectp2_y = 0.01;
    
    rectp3_x = -0.03;
    rectp3_y = -0.01;
    
    rectp4_x = 0.03;
    rectp4_y = -0.01;
    
    t = linspace(0,2*pi);
    
    rot_mat = [cos(orientation(i)) -sin(orientation(i)); sin(orientation(i)) cos(orientation(i))];
    
    xyoctagonpoints = rot_mat*[p1_x p2_x p3_x p4_x p5_x p6_x p7_x p8_x; p1_y p2_y p3_y p4_y p5_y p6_y p7_y p8_y];
    xycirclepoints = [cos(t);sin(t)];
    xyrectanglepoints = rot_mat*[rectp1_x rectp2_x rectp3_x rectp4_x; rectp1_y rectp2_y rectp3_y rectp4_y];

    barrelPatch(i).pointsX = [xyoctagonpoints(1,:)]' + xylocation(i,1);
    barrelPatch(i).pointsY = [xyoctagonpoints(2,:)]' + xylocation(i,2);
    barrelPatch(i+size(xylocation,1)).pointsX = [(r-.01)*xycirclepoints(1,:)]' + xylocation(i,1);
    barrelPatch(i+size(xylocation,1)).pointsY = [(r-.01)*xycirclepoints(2,:)]' + xylocation(i,2);
    barrelPatch(i+2*size(xylocation,1)).pointsX = [(r-.025)*xycirclepoints(1,:)]' + xylocation(i,1);
    barrelPatch(i+2*size(xylocation,1)).pointsY = [(r-.025)*xycirclepoints(2,:)]' + xylocation(i,2);
    barrelPatch(i+3*size(xylocation,1)).pointsX = [(r-.04)*xycirclepoints(1,:)]' + xylocation(i,1);
    barrelPatch(i+3*size(xylocation,1)).pointsY = [(r-.04)*xycirclepoints(2,:)]' + xylocation(i,2);
    barrelPatch(i+4*size(xylocation,1)).pointsX = [(r-.055)*xycirclepoints(1,:)]' + xylocation(i,1);
    barrelPatch(i+4*size(xylocation,1)).pointsY = [(r-.055)*xycirclepoints(2,:)]' + xylocation(i,2);    
    barrelPatch(i+5*size(xylocation,1)).pointsX = [xyrectanglepoints(1,:)]' + xylocation(i,1);
    barrelPatch(i+5*size(xylocation,1)).pointsY = [xyrectanglepoints(2,:)]' + xylocation(i,2);
end

