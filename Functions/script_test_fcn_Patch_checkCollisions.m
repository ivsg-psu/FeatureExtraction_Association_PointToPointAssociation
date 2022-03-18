% Script to generate vehicle path and collision data using a constant
% radius circular model of the vehicle behavior. Assumes the vehicle CG
% travels on the circle, defined as 

% Notation: px denotes points. For example, p0 is initial point, pc is the
% center point, etc.

% Set the getNewData flag to 1 to set new obstacles by mouse clicking
flag_getNewObst = 0;

% Try to add the path to the geometry library containing the function to
% find intersections between circles and line segments
addpath('~/Documents/MATLAB/PSU_GeometryLib/Functions/')
% Confirm that the function path was added successfully, or throw an error
if ~exist('fcn_geometry_findIntersectionLineSegmentWithCircle')
  error('Cannot find needed dependency fcn_geometry_findIntersectionLineSegmentWithCircle.m')
end

% Clear existing data if new data is desired
if 1 == flag_getNewObst
  clearvars -except flag*
end

% Enumerate the body positions and directions
LF = 1; RF = 2; RR = 3; LR = 4; X = 1; Y = 2;

% Vehicle trajectory information
vx = 20;                      % longitudinal speed (m/s)
R = 15;                       % path radius (m) with sign (+ left, - right)
p0 = [-8,2];                  % initial position of vehicle (m,m)
CoG0 = 4*pi*rand(1)-2*pi;     % initial course-over-ground of vehicle (rad)
a0 = (30*rand(1)-15)*pi/180;  % vehicle body slip angle (rad)
tf = 1.8*pi*abs(R)/vx;        % time horizon to plot (arbitrary)
% Create a vector of the trajectory information
x0 = [p0'; CoG0; a0; vx; R];
% Vehicle dimensional information
vehicle.dr = 8.2;       % CG-front bumper distance (m)
vehicle.df = 2.2;        % CG-rear bumper distance (m)
vehicle.w = 2.0;        % vehicle width (m)

% Reset ranges on angles to avoid having to handle multiple wraps later in
% the code
while CoG0 < 0
  CoG0 = CoG0 + 2*pi;
end
while CoG0 > 2*pi
  CoG0 = CoG0 - 2*pi;
end
while a0 < 0
  a0 = a0 + 2*pi;
end
while a0 > 2*pi
  a0 = a0 - 2*pi;
end

% Set up a figure for plotting the collision geometry
figure(1)
clf
hold on
grid on
axis equal

% Plot the vehicle CG at the initial location
plot(p0(1),p0(2),'k*')
% Plot a plus sign over the circle to make something akin to a CG symbol
plot(p0(1),p0(2),'k+','markersize',8)

% Calculate center point of vehicle trajectory circle
pc(1) = p0(1) + R*cos(CoG0+pi/2);
pc(2) = p0(2) + R*sin(CoG0+pi/2);

% Generate a circular trajectory for the given time horizon
del_t = 0.01;
t = (0:del_t:tf)';
theta = vx*t/R + CoG0 - sign(R)*pi/2;
N = length(t);
pv = zeros(N,2);
pv(:,1) = abs(R)*cos(theta) + pc(1);
pv(:,2) = abs(R)*sin(theta) + pc(2);
% Plot the circular vehicle CG trajectory
plot(pv(:,1),pv(:,2),'k-.')

% Calculate the various pertinent radii and corner points of the vehicle
[radii,vehicleBB,radiiFlags] = fcn_Patch_CalcCircularTrajectoryGeometry(x0,vehicle);

% Parse out which radii are which in the geometry
RinsideFront = sign(R)*radii(1);
RoutsideFront = sign(R)*radii(2);
Rmin = sign(R)*radii(6);
Rmax = sign(R)*radii(7);

% Output some descriptive information about the circular trajectory and the
% min and max radii cases
fprintf(1,"Min radius is %0.2f, max radius is %0.2f\n",Rmin,Rmax);
locations = {'LF','RF','RR','LR','tangent'};
fprintf(1,"Min radius on the %s point, max radius on the %s point\n",...
  locations{radiiFlags(1)},locations{radiiFlags(2)});

% Calculate the offset to the angular position based on the back end of the
% vehicle (in order to plot the extra portion of the clearance curves)
rear_angles = atan2(vehicleBB(RR:LR,Y)-pc(Y),vehicleBB(RR:LR,X)-pc(X));
% Depending on the travel direction, determine the angle to the rearmost
% point on the vehicle bounding obx
if R > 0
  [~,start_idx] = min(rear_angles);
  start_angle = rear_angles(start_idx);
  % Correct for cases where wrapping has caused a mis-ordering of the start
  % and end angles
  if theta(1) - start_angle > 2*pi
    start_angle = start_angle + 2*pi;
  end
  % Create the vector of points to plot the various radii
  theta_extra = start_angle:pi/180:CoG0 - pi/2;
else
  [~,start_idx] = max(rear_angles);
  start_angle = rear_angles(start_idx);
  % Correct for cases where wrapping has caused a mis-ordering of the start
  % and end angles
  if theta(1) - start_angle > pi
    start_angle = start_angle + 2*pi;
  end
  % Create the vector of points to plot the various radii
  theta_extra = start_angle:-pi/180:CoG0 + pi/2;
end
% Plot the min and max radii along with the radii for the front corners of
% the vehicle
plot(abs(RinsideFront)*cos(theta_extra) + pc(1),abs(RinsideFront)*sin(theta_extra) + pc(2),'-.','color',[0.7 0.7 0.7]);
plot(abs(RinsideFront)*cos(theta) + pc(1),abs(RinsideFront)*sin(theta) + pc(2),'-','color',[0.7 0.7 0.7]);
plot(abs(RoutsideFront)*cos(theta_extra) + pc(1),abs(RoutsideFront)*sin(theta_extra) + pc(2),'-.','color',[0.7 0.7 0.7]);
plot(abs(RoutsideFront)*cos(theta) + pc(1),abs(RoutsideFront)*sin(theta) + pc(2),'-','color',[0.7 0.7 0.7]);
plot(abs(Rmin)*cos(theta_extra) + pc(1),abs(Rmin)*sin(theta_extra) + pc(2),'r-.');
plot(abs(Rmin)*cos(theta) + pc(1),abs(Rmin)*sin(theta) + pc(2),'r-');
plot(abs(Rmax)*cos(theta_extra) + pc(1),abs(Rmax)*sin(theta_extra) + pc(2),'b-.');
plot(abs(Rmax)*cos(theta) + pc(1),abs(Rmax)*sin(theta) + pc(2),'b-');

% Plot the initial vehicle location
plotVehicleBB(p0,CoG0+a0,vehicle,1)

% Create some obstacles in the global coordinates
if 1 == flag_getNewObst
  %pause(); % Used for zooming the plot to create more accurate obstacles, if desired
  Nobstacles = 1;
  % Create an template obstacle structure array to fill
  obstacles = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});
  
  for obstacleInd = 1:Nobstacles
    % Add an obstacle by clicking to generate the points and then
    % turning the points into a patch object
    [xobst,yobst] = ginput;
    obstacles(obstacleInd).pointsX = xobst;
    obstacles(obstacleInd).pointsY = yobst;
    obstacles(obstacleInd).color = [0.4 0.4 0.4];
  end
end

% Plot the patch object over the trajectory lines using the patch plotting
% utility
[h,hpts] = fcn_Patch_plotPatch(obstacles,1);
% Update the axis limits to see the complete scenario
axis auto

%% Check for a collision at the initial configuration

% Set a flag that will be set if any of the obstacles yield a collision
initial_collision_flag = 0;

% Loop through each obstacle, checking for collisions with the vehicle
for obstacleInd = 1:Nobstacles
  % Use the polyxpoly function to check for overlaps between the car
  % bounding box and obstacle polyshapes
  [XI,YI] = polyxpoly(vehicleBB([1:end 1],X),vehicleBB([1:end 1],Y),obstacles(obstacleInd).pointsX([1:end 1]),obstacles(obstacleInd).pointsY([1:end 1]));
  
  % If there are entries in the intersection vectors, there is a
  % collision. Checking XI only since XI and YI are the same size.
  if ~isempty(XI)
    % Since the intersection vectors were non-empty, there is a
    % collision. Set the flag appropriately.
    initial_collision_flag = 1;
    % Output the collision notification
    fprintf(1,'Vehicle initial position is in conflict with obstacle %d.\n',obstacleInd);
    % Determine the overlap between the vehicle and obstacle
    collisionPoly = intersect(polyshape(vehicleBB([1:end 1],X),vehicleBB([1:end 1],Y)),polyshape(obstacles(obstacleInd).pointsX([1:end 1]),obstacles(obstacleInd).pointsY([1:end 1])));
    % Plot the overlapped region in red for the user's review
    plot(collisionPoly,'facecolor','red')
  end
end

%% Determine the nearest collision

% Run this check only if the vehicle does not begin in a state of collision
if 0 == initial_collision_flag
  tic     % Start timer
  % Check for straight-line motion
  if Inf == abs(R)
    error('This script does not handle straight line motion. Use script_test_fcn_Patch_checkStraightCollisions instead.');
  else
    % Run the main collision detection algorithm
    [collFlags,collTime,collAngle,collLoc,clearance,bodyCollLoc] = fcn_Patch_checkCollisions(x0,vehicle,obstacles)
  end
  ET = toc;   % Stop timer and store the elapsed time
  
  fprintf(1,'Determined collision geometry for %d objects in %0.3f seconds.\n',length(collFlags),ET);
  
  % Plot the car body
  Rabs = abs(R);
  frontEdgeOffset = sign(R)*pi/2;
  
  % Iterate over all of objects
  for obstacleInd = 1:Nobstacles
    % Determine the vehicle heading at the point of collision or near miss
    % using the angle of the initial vehicle position, theta(1)
    collHeading = theta(1) + sign(R)*collAngle(obstacleInd) + (pi/2+a0);
    % Determine the CG location at the point of collision or near miss
    collCG = [R*cos(theta(1) + sign(R)*collAngle(obstacleInd)) ...
              R*sin(theta(1) + sign(R)*collAngle(obstacleInd))] + pc;
    % Plot the vehicle body in the position of the collision or near miss
    plotVehicleBB(collCG,collHeading,vehicle,1)
    
    % Plot the collision or closest clearance point
    plot(collLoc(obstacleInd,1),collLoc(obstacleInd,2),'r*')
    
    % Report the findings, either for a collision or near miss
    if 1 == collFlags(obstacleInd)
      fprintf(1,'  For object %d, found a collision with TTC of %0.2f seconds at (%0.2f,%0.2f)\n',obstacleInd,collTime(obstacleInd),collLoc(obstacleInd,1),collLoc(obstacleInd,2));
    else
      fprintf(1,'  For object %d, no collision detected. Smallest clearance distance is %0.2f units at (%0.2f,%0.2f), occurring at %0.2f seconds.\n',obstacleInd,clearance(obstacleInd),collLoc(obstacleInd,1),collLoc(obstacleInd,2),collTime(obstacleInd));
    end
  end
  
  [~,firstColl] = conditionalMin(collTime,collFlags,'== 1');
  if isempty(firstColl)
    firstColl = 1;
  end
  % Activate figure(1) in case it was de-selected
  figure(1)
  % Grab the handle to the axis object of the current figure
  a1 = gca;
  % Create or address figure 2
  f2 = figure(2);
  % Clear the contents of figure 2, if any
  clf
  % Copy all of the contents of the axis object in figure 1 to figure 2
  a2 = copyobj(a1,f2);
  % Address figure 2
  figure(2)
  % Set the axes fairly tightly around the collision location
  axis([collLoc(firstColl,1)-2 collLoc(firstColl,1)+2 collLoc(firstColl,2)-2 collLoc(firstColl,2) + 2]);
end

function plotVehicleBB(pCG,heading,vehicle,figHandle)

% Enumerate the body positions and directions
LF = 1; RF = 2; RR = 3; LR = 4; X = 1; Y = 2;

% Address the specified figure
figure(figHandle)
% Plot the CG location
plot(pCG(X),pCG(Y),'ko')
plot(pCG(X),pCG(Y),'k+')

% Compute the locations of the bounding box points by adding vectors to the
% CG. To move fore/aft, add df/dr at the heading angle or the heading angle
% - pi. To move left/right, add w/2 at the heading angle plus/minus pi/2.
vehicleBB(LF,X) = pCG(X) + vehicle.df*cos(heading)    + vehicle.w/2*cos(heading + pi/2);
vehicleBB(LF,Y) = pCG(Y) + vehicle.df*sin(heading)    + vehicle.w/2*sin(heading + pi/2);
vehicleBB(RF,X) = pCG(X) + vehicle.df*cos(heading)    + vehicle.w/2*cos(heading - pi/2);
vehicleBB(RF,Y) = pCG(Y) + vehicle.df*sin(heading)    + vehicle.w/2*sin(heading - pi/2);
vehicleBB(RR,X) = pCG(X) + vehicle.dr*cos(heading-pi) + vehicle.w/2*cos(heading - pi/2);
vehicleBB(RR,Y) = pCG(Y) + vehicle.dr*sin(heading-pi) + vehicle.w/2*sin(heading - pi/2);
vehicleBB(LR,X) = pCG(X) + vehicle.dr*cos(heading-pi) + vehicle.w/2*cos(heading + pi/2);
vehicleBB(LR,Y) = pCG(Y) + vehicle.dr*sin(heading-pi) + vehicle.w/2*sin(heading + pi/2);
% Plot the vehicle bounding box
plot(vehicleBB([1:end 1],X),vehicleBB([1:end 1],Y),'b-','linewidth',1);
end
