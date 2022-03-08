function [trajectoryRadii,boundingPoints,radiiFlags] = fcn_Patch_CalcCircularTrajectoryGeometry(x0,vehicle)
% fcn_Patch_CalcCircularTrajectoryGeometry
% Evaluates the initial position and orientation of a vehicle that will
% traverse a circular trajectory. Returns relevant radii for the trajectory
% as well as the initial body bounding box points.
%
% ASSUMPTIONS:
%       1) The vehicle moves at constant speed along a circular trajectory.
%       2) The vehicle is represented by a bounding rectangle.
%       3) The vehicle slip angle is limited to [-90,90] degrees
%       4) A positive trajectory radius indicates CCW travel. Negative
%       indicates CW.
%
% FORMAT:
%
%       [time,angle,location,clearance,bodyLoc] = fcn_Patch_checkCollisions(x0,vehicle,patchArray)
%
% INPUTS:
%
%      x0: a 6 x 1 vector containing the starting (x,y) coordinates of the
%           vehicle, the initial course over ground (CoG), the body slip
%           angle, the longitudinal vehicle speed, and the signed
%           trajectory radius in (m,m), radians, radians, m/s, and m.
%           (Speed is not needed here, but is included for consistency with
%           other vehicle trajectory functions.)
%      vehicle: a structure containing the vehicle properties, which must
%           include fields df, dr, w for the vehicle CG-front bumper
%           distance, CG-rear bumper distance, and body width, in meters,
%           respectively.
%
% OUTPUTS:
%
%      trajectoryRadii: a 5 x 1 vector representing the relevant radii for
%           a vehicle traveling around a circular trajectory, with the
%           following elements:
%           (2) RinsideFront: a scalar value representing the radius upon
%               which the inside front wheel travels (sometimes coincident
%               with Rmin, depending on the vehicle slip angle)
%           (3) RinsideRear: a scalar value representing the radius upon
%               which the inside rear wheel travels (sometimes coincident
%               with Rmin, depending on the vehicle slip angle)
%           (4) RoutsideFront: a scalar value representing the radius upon
%               which the outside front wheel travels (sometimes coincident
%               with Rmax, depending on the vehicle slip angle)
%           (5) RoutsideRear: a scalar value representing the radius upon
%               which the outside rear wheel travels (sometimes coincident
%               with Rmax, depending on the vehicle slip angle)
%           (6) RinsideEdge: a scalar value representing the radius that is
%               tangent to the inside edge of the vehicle and is sometimes
%               the minimum clearance radius
%           (6) Rmin: smallest radius on the circular trajectory that the 
%               vehicle will "clear" and not collide with an obstacle
%               placed inside this radius
%           (7) Rmax: largest radius on the circular trajectory that the 
%               vehicle will "clear" and not collide with an obstacle
%               placed outside of this radius
%      boundingPoints: a 4 x 2 vector representing the locations of the
%           corners of the bounding box for the vehicle, in order of LF,
%           LR, RF, RR
%
% DEPENDENCIES:
%
%       No dependencies.
%
% EXAMPLES:
%
%       See the script: script_test_fcn_Patch_CalcCircularTrajectoryGeometry.m for a full test
%       suite.
%
% This function was written by C. Beal
% Questions or comments? cbeal@bucknell.edu

% Revision history:
%     2022_03_08
%     -- wrote the code

flag_do_debug = 1; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    fh = gcf;   % grab the handle for the current figure
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
    if nargin ~= 2
        error('Incorrect number of input arguments')
    end
    
    % Check the trajectory input
    if size(x0,1) ~= 6
        error('Vehicle trajectory missing elements, cannot calculate collisions.')
    end
    
    % Make sure the trajectory isn't degenerate with a radius of zero
    if 0 == x0(6)
        error('Vehicle trajectory radius must have magnitude greater than zero');
    end
    
    % Make sure the vehicle slip angle is not being 90 degrees in magnitude
    if pi/2 < abs(x0(4))
        error('Vehicle slip angle must have magnitude less than pi/2 (90 degrees)');
    end    
    
    % Check the vehicle structure input to make sure that the dimensions a,
    % b, and d are all supplied
    if ~all(isfield(vehicle,{'df','dr','w'}))
        error('One or more necessary vehicle dimensions missing. Check inputs.')
    end
    if isfield(vehicle,'df') && isempty(vehicle.df)
        error('CG-front axle distance empty.')
    end
    if isfield(vehicle,'dr') && isempty(vehicle.dr)
        error('CG-rear axle distance empty.')
    end
    if isfield(vehicle,'w') && isempty(vehicle.w)
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
CoG0 = x0(3);   % Initial course over ground (CoG) of the vehicle
a0 = x0(4);     % Vehicle body slip angle
v0 = x0(5);     % Initial (and constant) vehicle speed
R = x0(6);      % Magnitude of the radius of the circular trajectory

% Enumerate the body bounding box points for easier reading of the code
LF = 1;
RF = 2;
RR = 3;
LR = 4;
% Enumerate the directions in the matrices for easier reading of the code
X = 1;
Y = 2;

% Calculate center point of vehicle trajectory circle
pc(1) = p0(X) + R*cos(CoG0+pi/2);
pc(2) = p0(Y) + R*sin(CoG0+pi/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First, calculate the locations of each of the corner points of the
% vehicle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_do_debug
figure(2)
clf
hold on
grid on
axis equal
end

% Set up the vehicle corner points in standard position (facing along the
% positive x-axis)
boundingPointsNorm = [vehicle.df vehicle.w/2; vehicle.df -vehicle.w/2;...
    -vehicle.dr -vehicle.w/2; -vehicle.dr vehicle.w/2];
% Plot the bounding box in standard position and orientation
if flag_do_debug
plot(boundingPointsNorm([1:end 1],X),boundingPointsNorm([1:end 1],Y),'k-o');
end

% Now, create a matrix to rotate points by the CoG and slip angle
rotMat = [cos(CoG0+a0) sin(CoG0+a0); -sin(CoG0+a0) cos(CoG0+a0)];
% Multiply each X,Y bounding box corner by the rotation matrix
boundingPointsRotated = (rotMat'*boundingPointsNorm')';
% Plot the new bounding box
if flag_do_debug
plot(boundingPointsRotated([1:end 1],X),boundingPointsRotated([1:end 1],Y),'b-o');
end

% Shift the bounding points so that the vehicle is centered at the initial
% point
boundingPoints(:,X) = boundingPointsRotated(:,X) + p0(X);
boundingPoints(:,Y) = boundingPointsRotated(:,Y) + p0(Y);
% Plot the new bounding box
if flag_do_debug
plot(boundingPoints([1:end 1],X),boundingPoints([1:end 1],Y),'r-o');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next, calculate the radii associated with each location on the vehicle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the radius associated with each of the corners of the vehicle
% bounding box
trajectoryRadii(LF) = norm(boundingPoints(LF,:)-pc);
trajectoryRadii(RF) = norm(boundingPoints(RF,:)-pc);
trajectoryRadii(RR) = norm(boundingPoints(RR,:)-pc);
trajectoryRadii(LR) = norm(boundingPoints(LR,:)-pc);

% Compute the alpha associated with the min radius
if R >= 0
    pa = boundingPoints(LR,:); pb = boundingPoints(LF,:);
else
    pa = boundingPoints(RR,:); pb = boundingPoints(RF,:);
end

% Find the fraction of the distance from pa to pb at which the minimum
% radius is tangent to the line segment a-b
alpha = -((pa(X)-pc(X))*(pb(X)-pa(X)) + (pa(Y)-pc(Y))*(pb(Y)-pa(Y)))/((pb(X)-pa(X))^2 + (pb(Y)-pa(Y))^2);
% The fraction isn't bounded and can be greater than 1 or negative, which
% would indicate an intersection outside of the line segment limits, so
% only handle this case if it happens along the actual side of the vehicle.
if alpha < 1 && alpha > 0
    % Store the tangent point between the inner edge and the minimum radius
    boundingPoints(5,:) = pa + alpha*(pb-pa);
    % Add the point to the plot, if desired
    if flag_do_debug
        plot(boundingPoints(5,X),boundingPoints(5,Y),'r*')
    end
    % Calculate and store the radius associated with this tangent point
    trajectoryRadii(5) = norm(boundingPoints(5,:)-pc);
else
    % If the tangent point lies outside of the vehicle, store a NaN for
    % this point
    trajectoryRadii(5) = NaN;
end

% Calculate the minimum radius, including the potential tangent point. Need
% to use nanmin because of the possibility that the tangent point does not
% exist. Store the point that yielded the minimum radius for potential
% logical processing outside of this function.
[trajectoryRadii(6),radiiFlags(1)] = nanmin(trajectoryRadii(1:5));
% Calculate the maximum radius and store the point that yielded the maximum
% radius for potential logical processing outside of this function.
[trajectoryRadii(7),radiiFlags(2)] = max(trajectoryRadii(1:4));

if flag_do_debug
    % If desired, output the min/max radii and locations they were found
    fprintf(1,"Min radius is %0.2f, max radius is %0.2f\n",trajectoryRadii(6),trajectoryRadii(7));
    locations = {'LF','RF','RR','LR','tangent'};
    fprintf(1,"Min radius on the %s point, max radius on the %s point\n",locations{radiiFlags(1)},locations{radiiFlags(2)});
    
    % Set the focus back to the original figure to avoid screwing up
    % calling scripts/functions. This needs to come after any other
    % plotting commands in this function to avoid shifting back or plotting
    % extra things in the user's original figure.
    figure(fh);
end
