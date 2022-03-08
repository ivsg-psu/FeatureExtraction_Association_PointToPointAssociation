function [collFlag,time,location,clearance,bodyLoc] = fcn_Patch_checkStraightCollisions(x0,vehicle,patchArray)
% fcn_Patch_checkStraightCollisions
% Evaluates a straight vehicle trajectory against a series of patches to
% determine whether there will be collisions between the vehicle and the
% outline of any of the patch objects.
%
% ASSUMPTIONS:
%       1) The vehicle moves at constant speed along a straight trajectory.
%       2) The vehicle travels with zero sideslip angle.
%       3) The vehicle is represented by a bounding rectangle.
%
% FORMAT:
%
%       [flags,time,location,clearance,bodyLoc] = fcn_Patch_checkStraightCollisions(x0,vehicle,patchArray)
%
% INPUTS:
%
%      x0: a 4 x 1 vector containing the starting (x,y) coordinates of the
%           vehicle, the initial heading, and the longitudinal vehicle
%           speed, in (m,m), radians, m/s, and m.
%      vehicle: a structure containing the vehicle properties, which must
%           include fields df, dr, w for the vehicle CG-front bumper
%           distance, CG-rear bumper distance, and body width, in meters,
%           respectively.
%      patchArray: a structure array defining the objects with which the
%           vehicle could potentially collide
%
% OUTPUTS:
%
%      collFlag: an N x 1 vector of flags denoting whether there is a
%           collision with each of the objects
%      time: an N x 1 vector of collision times, where N is the number of
%           patch objects. Elements of the time vector will be set to Inf
%           if there is no overlap with the vehicle path.
%      location: a N x 2 vector of collision locations, where N is the
%           number of patch objects and the columns are the x and y
%           coordinates of the collision. Elements of the location matrix
%           will be set to NaN if there is no overlap with the vehicle
%           path.
%       clearance: an N x 1 vector of minimum clearance distances between
%           the vehicle and the patch objects, where N is the number of
%           patch objects. Elements of the clearance vector will be set to
%           NaN if there is a collision with the object.
%       bodyLoc: an N x 2 vector of collision locations in vehicle body
%           fixed coordinates, where N is the number of patch objects and
%           the columns are the x and y coordinates of the collision.
%           Elements of the location matrix will be set to NaN if there is
%           no overlap with the vehicle path
%
% DEPENDENCIES:
%
%      ## NOT CURRENTLY USED: fcn_Patch_checkInputsToFunctions
%
% EXAMPLES:
%
%       See the script: script_test_fcn_Patch_checkCollisions.m for a full test
%       suite.
%
% This function was written by C. Beal
% Questions or comments? cbeal@bucknell.edu

% Revision history:
%     2022_03_03
%     -- adapted the code from the code for circular trajectories

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
    if nargin ~= 3
        error('Incorrect number of input arguments')
    end
    
    % Check to see that there were any points provided
    if isempty(patchArray)
        warning('Empty array of objects, nothing to do.');
        return
    end
    
    % Check the trajectory input
    if size(x0,1) ~= 4
        error('Vehicle trajectory missing elements, cannot calculate collisions.')
    end
    
    % Check the vehicle structure input to make sure that the dimensions a,
    % b, and d are all supplied
    if ~all(isfield(vehicle,{'df','dr','w'}))
        error('One or more necessary vehicle dimensions missing. Check inputs.')
    end
    if isfield(vehicle,'df') && isempty(vehicle.df)
        error('CG-front axle distance empty.')
    end
    if isfield(vehicle,'dr') && isempty(vehicle.df)
        error('CG-rear axle distance empty.')
    end
    if isfield(vehicle,'w') && isempty(vehicle.df)
        error('Vehicle width empty.')
    end
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

% Break out some variables for easier referencing
p0 = x0(1:2);   % Initial location of the vehicle
h0 = x0(3);     % Initial heading of the vehicle
v0 = x0(4);     % Initial (and constant) vehicle speed

% Determine the number of patches to check
Npatches = length(patchArray);

% Define the rotation matrix that rotates the vehicle to point to the
% positive x-axis in order to rotate the obstacles by the same angle
rotMat = [cos(h0) sin(h0); -sin(h0) cos(h0)];

% Iterate over all of the patches to check
for patchInd = 1:Npatches
    % Grab the X,Y points out of the patch structure
    xobst = patchArray(patchInd).pointsX;
    yobst = patchArray(patchInd).pointsY;
    
    % Shift the X,Y points by the vehicle offset
    xobst = xobst - p0(1);
    yobst = yobst - p0(2);
    
    % Rotate the scenario/points around the origin to align the vehicle
    % travel with the x-axis to turn this into a "standard" problem
    rotObst = (rotMat*[xobst yobst]')';
    
    % Determine any obstacles that lie within or across the vehicle path on
    % the positive the x-axis between the back of the vehicle (-vehicle.dr)
    % and some large distance (1e6)
    [xColl,yColl] = polyxpoly([-vehicle.dr 1e6 1e6 -vehicle.dr -vehicle.dr],...
        [-1 -1 1 1 -1]*vehicle.w/2,rotObst([1:end 1],1),rotObst([1:end 1],2));
    
    % Determine any vertices that fall within the trajectory that are
    % within the width of the vehicle on the x-axis between the back of the
    % vehicle (-vehicle.dr) and some large distance (1e6)
    collVertexInds = find(abs(rotObst(:,2)) < 1 & rotObst(:,1) > -vehicle.dr);
    
    % Add any vertices that are within the path to the list of collision
    % locations that will be evaluated to find the earliest one on the path
    xColl = [xColl; rotObst(collVertexInds,1)];
    yColl = [yColl; rotObst(collVertexInds,2)];
    
    % Determine whether there is a collision at all
    if ~isempty(xColl)
        % There is a collision if we've gotten here, so set the flag to
        % indicate a collision and set the clearance to NaN
        collFlag(patchInd) = 1;
        clearance(patchInd) = NaN;
        
        % If there is a collision immediately, find the farthest one
        % toward the back of the vehicle (that would have happened
        % earliest)
        [~,closestCollInd] = min(xColl);
        % Check to see if the collision happens immediately (already in
        % collision state)
        if xColl(closestCollInd) < vehicle.df
            % Fill the output variable for the collision location, rotated
            % and translated back to the original coordinate system
            location(patchInd,:) = ([xColl(closestCollInd) yColl(closestCollInd)]*rotMat)' + p0;
            % Call out the location on the vehicle body as the point on the
            % width of the front of the vehicle that matches the obstacle
            bodyLoc(patchInd,:) = [vehicle.df yColl(closestCollInd)];
            % Compute the time of the collision (negative since it already
            % happened)
            time(patchInd) = (xColl(closestCollInd)-vehicle.df)/v0;
        else
            % Determine the time to reach the collision location that is in
            % front of the vehicle
            time(patchInd) = (xColl(closestCollInd)-vehicle.df)/v0;
            % The object will collide with the front of the vehicle, but
            % use the y-coordinate of the collision to determine where on
            % the front of the vehicle
            bodyLoc(patchInd,:) = [vehicle.df yColl(closestCollInd)];
            % Fill the output variable for the collision location, rotated
            % and translated back to the original coordinate system
            location(patchInd,:) = ([xColl(closestCollInd) yColl(closestCollInd)]*rotMat)' + p0;
        end
    else
        % Since there is not a collision, set the flag appropriately and
        % set the body collision location to NaN
        collFlag(patchInd) = 0;
        % Now find the minimum clearance of all of the points on the
        % obstacle
        [clearance(patchInd),minClearanceInd] = min(abs(rotObst(:,2))-vehicle.w/2);
        % Determine whether the near miss is on the right or left front
        % corner of the vehicle and store the result
        bodyLoc(patchInd,:) = [vehicle.df vehicle.w/2*sign(rotObst(minClearanceInd,2))];
        % Determine the time to get to the minimum clearance location with
        % the front of the vehicle
        time(patchInd) = (rotObst(minClearanceInd,1)-vehicle.df)/v0;
        % Store the location of the minimum clearance location, rotated
        % and translated back to the original coordinate system
        location(patchInd,:) = (rotObst(minClearanceInd,:)*rotMat)' + p0;
    end  
end