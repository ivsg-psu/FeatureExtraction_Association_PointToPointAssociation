% Script to generate vehicle path and collision data using a constant
% radius circular model of the vehicle behavior. Currently assumes the
% vehicle CG travels on the circle and that there is zero sideslip angle.
% Could be extended to include vehicle sideslip, though the approach
% would have to use either max sideslip angle or change the bounds of the
% vehicle as it goes around.

% Notation: px denotes points. For example, p0 is initial point, pc is the
% center point, etc.

% Set the getNewData flag to 1 to set new obstacles by mouse clicking
flag_getNewData = 1;
flag_loadData = 0;

addpath('~/Documents/MATLAB/PSU_GeometryLib/Functions/')

% Clear existing data if new data is desired
if 1 == flag_getNewData || 1 == flag_loadData
    clearvars -except flag*
%     % Reset the flag for later use since clearvars will clear it
%     flag_getNewData = 1;
    if 1 == flag_loadData
        load test_obstacle
        flag_getNewData = 0;
    end
end

% Enumerate the body positions and directions
LF = 1; RF = 2; RR = 3; LR = 4; X = 1; Y = 2;

% Vehicle trajectory information
vx = 20;        % longitudinal speed (m/s)
R = 15;         % path radius (m) with sign (+ left, - right)
p0 = [10,-5];     % initial position of vehicle (m,m)
h0 = pi/2;     % initial heading of vehicle (rad)
a0 = 15*pi/180; % vehicle body slip angle (rad)
tf = 1.9*pi*abs(R)/vx;         % time horizon to check (s)
% Create a vector of the trajectory information
x0 = [p0'; h0; a0; vx; R];
% Vehicle dimensional information
vehicle.dr = 2.2;       % CG-front bumper distance (m)
vehicle.df = 8;        % CG-rear bumper distance (m)
vehicle.w = 2.0;        % vehicle width (m)

% Plot the points
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
pc(1) = p0(1) + R*cos(h0+pi/2);
pc(2) = p0(2) + R*sin(h0+pi/2);

% Generate a trajectory for the given time horizon
del_t = 0.01;
t = (0:del_t:tf)';
theta = vx*t/R + h0 - pi/2;
N = length(t);
pv = zeros(N,2);
pv(:,1) = R*cos(theta) + pc(1);
pv(:,2) = R*sin(theta) + pc(2);
% And plot the trajectory
plot(pv(:,1),pv(:,2),'k-.')

% Calculate the various pertinent radii and corner points of the vehicle
[radii,vehicleBB,radiiFlags] = fcn_Patch_CalcCircularTrajectoryGeometry(x0,vehicle);

% Parse out which radii are which
Rinside = sign(R)*radii(1);
RoutsideFront = sign(R)*radii(2);
RoutsideRear = sign(R)*radii(3);
Rmin = sign(R)*radii(6);
Rmax = sign(R)*radii(7);

% The radii flags can be used to determine which edges of the vehicle are
% of concern for calculating collisions. Possibilities are left, front, and
% right (radiiFlags = [5,3]/[tangent,RR] or [4,3]/[LR,RR]), left and front
% (radiiFlags = [5,2]/[tangent,RF] or [4,2]/[LR,RF]), or front and right
% (radiiflags = [1,3]/[LF,RR])
fprintf(1,"Min radius is %0.2f, max radius is %0.2f\n",Rmin,Rmax);
locations = {'LF','RF','RR','LR','tangent'};
fprintf(1,"Min radius on the %s point, max radius on the %s point\n",...
    locations{radiiFlags(1)},locations{radiiFlags(2)});

% Calculate the offset to the angular position based on the back end of the
% vehicle (in order to plot the extra portion of the clearance curves)
rear_offset = 1.1*sign(R)*atan2(-vehicle.dr,sign(R)*R+vehicle.w/2);
theta_arcs = linspace(theta(1)+rear_offset,theta(end),100);
plot(Rinside*cos(theta_arcs) + pc(1),Rinside*sin(theta_arcs) + pc(2),'-','color',[0.7 0.7 0.7]);
plot(RoutsideFront*cos(theta_arcs) + pc(1),RoutsideFront*sin(theta_arcs) + pc(2),'-','color',[0.7 0.7 0.7]);
plot(Rmin*cos(theta_arcs) + pc(1),Rmin*sin(theta_arcs) + pc(2),'r-');
plot(Rmax*cos(theta_arcs) + pc(1),Rmax*sin(theta_arcs) + pc(2),'b-');

% Plot the initial vehicle location
plotVehicleBB(p0,h0+a0,vehicle,1)

% Create some obstacles in the global coordinates
if 1 == flag_getNewData
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
        [collFlags,collTime,collAngle,collLoc,clearance,bodyCollLoc] = fcn_Patch_checkCollisions(x0,vehicle,obstacles);
    end
    ET = toc;   % Stop timer and store the elapsed time
    
    fprintf(1,'Determined collision geometry for %d objects in %0.3f seconds.\n',length(collFlags),ET);
    
    % Plot the car body
    Rabs = abs(R);
    frontEdgeOffset = sign(R)*pi/2;
  %%  
    for collInd = 1:length(collFlags)
        % Plot the vehicle body in the position of the collision or near
        % miss
        collHeading = theta(1) + sign(R)*collAngle(collInd) + (pi/2+a0);
        collCG = [R*cos(theta(1) + sign(R)*collAngle(collInd)) R*sin(theta(1) + sign(R)*collAngle(collInd))] + pc;
        plotVehicleBB(collCG,collHeading,vehicle,1)
        
        % Plot the collision or closest clearance point
        plot(collLoc(collInd,1),collLoc(collInd,2),'r*')
        
        if 1 == collFlags(collInd)
            fprintf(1,'  For object %d, found a collision with TTC of %0.2f seconds at (%0.2f,%0.2f)\n',collInd,collTime(collInd),collLoc(collInd,1),collLoc(collInd,2));
        else
            fprintf(1,'  For object %d, no collision detected. Smallest clearance distance is %0.2f units at (%0.2f,%0.2f), occurring at %0.2f seconds.\n',collInd,clearance(collInd),collLoc(collInd,1),collLoc(collInd,2),collTime(collInd));
        end
    end
    
    [~,firstColl] = conditionalMin(collTime,collFlags,'== 1');
    
%     % Create another copy of the figure
%     figure(1)
%     a1 = gca;
%     f2 = figure(2);
%     clf
%     a2 = copyobj(a1,f2);
%     figure(2)
%     axis([collLoc(firstColl,1)-2 collLoc(firstColl,1)+2 collLoc(firstColl,2)-2 collLoc(firstColl,2) + 2]);
    % Create a legend
    %legend('Vehicle Start Point','Vehicle Trajectory','Inner Vehicle Bound','Outer Vehicle Bound','Vehicle CG','Vehicle Outline at Start','Obstacle','Obstacle Vertices','Vehicle Outline at Collision','Collision Point','location','best')
end

function plotVehicleBB(pCG,heading,vehicle,figHandle)

% Enumerate the body positions and directions
LF = 1; RF = 2; RR = 3; LR = 4; X = 1; Y = 2;

% Address the specified figure
figure(figHandle)
% Plot the CG location
plot(pCG(X),pCG(Y),'ko')
plot(pCG(X),pCG(Y),'k+')

% Compute the locations of the bounding box points by adding vectors to the CG. To move fore/aft, add df/dr at the heading angle or the heading angle - pi. To move left/right, add w/2 at the heading angle plus/minus pi/2.
vehicleBB(LF,X) = pCG(X) + vehicle.df*cos(heading)    + vehicle.w/2*cos(heading + pi/2);
vehicleBB(LF,Y) = pCG(Y) + vehicle.df*sin(heading)    + vehicle.w/2*sin(heading + pi/2);
vehicleBB(RF,X) = pCG(X) + vehicle.df*cos(heading) + vehicle.w/2*cos(heading - pi/2);
vehicleBB(RF,Y) = pCG(Y) + vehicle.df*sin(heading) + vehicle.w/2*sin(heading - pi/2);
vehicleBB(RR,X) = pCG(X) + vehicle.dr*cos(heading-pi) + vehicle.w/2*cos(heading - pi/2);
vehicleBB(RR,Y) = pCG(Y) + vehicle.dr*sin(heading-pi) + vehicle.w/2*sin(heading - pi/2);
vehicleBB(LR,X) = pCG(X) + vehicle.dr*cos(heading-pi)    + vehicle.w/2*cos(heading + pi/2);
vehicleBB(LR,Y) = pCG(Y) + vehicle.dr*sin(heading-pi)    + vehicle.w/2*sin(heading + pi/2);
% Plot the vehicle bounding box
plot(vehicleBB([1:end 1],X),vehicleBB([1:end 1],Y),'b-','linewidth',1);
end
