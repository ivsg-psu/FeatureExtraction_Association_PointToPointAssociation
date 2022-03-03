% Script to generate vehicle path and collision data using a constant
% radius circular model of the vehicle behavior. Currently assumes the
% vehicle CG travels on the circle and that there is zero sideslip angle.
% Could be extended to include vehicle sideslip, though the approach
% would have to use either max sideslip angle or change the bounds of the
% vehicle as it goes around.

% Notation: px denotes points. For example, p0 is initial point, pc is the
% center point, etc.

getNewData = 1;
if 1 == getNewData
    clearvars
    getNewData = 1;
end

% Vehicle trajectory information
vx = 20;        % longitudinal speed (m/s)
R = 10;         % path radius (m) with sign (+ left, - right)
p0 = [0,0];     % initial position of vehicle (m,m)
h0 = -pi/2;      % initial heading of vehicle (rad)
tf = 2.36;         % time horizon to check (s)
% Create a vector of the trajectory information
x0 = [p0'; h0; vx; R];
% Vehicle dimensional information
vehicle.a = 1.8;       % CG-front axle distance (m)
vehicle.b = 2.2;        % CG-rear axle distance (m)
vehicle.d = 2.0;        % vehicle width (m)

% Plot the points
figure(1)
clf
hold on
grid on
axis equal
plot(p0(1),p0(2),'k*')

% Calculate center point of vehicle trajectory circle
pc(1) = p0(1) + R*cos(h0+pi/2);
pc(2) = p0(2) + R*sin(h0+pi/2);

% Generate a trajectory for the given time horizon
del_t = 0.01;
t = (0:del_t:tf)';
theta = vx*t/R;
N = length(t);
pv = zeros(N,2);
pv(:,1) = R*cos(theta+h0-pi/2) + pc(1);
pv(:,2) = R*sin(theta+h0-pi/2) + pc(2);
% And plot the trajectory
plot(pv(:,1),pv(:,2),'k-.')

% Determine the bounding radii for all portions of the vehicle
Rmin = R - vehicle.d/2*sign(R);
Rinner = sign(R)*sqrt((R-vehicle.d/2*sign(R))^2 + vehicle.a^2);
Router = sign(R)*sqrt((R+vehicle.d/2*sign(R))^2 + vehicle.a^2);
Rmax = sign(R)*sqrt((R+vehicle.d/2*sign(R))^2 + max(vehicle.a,vehicle.b)^2);
travel_offset = h0 - pi/2;
theta_max_offset = sign(R)*atan2(-vehicle.b,sign(R)*R+vehicle.d/2);
plot(Rinner*cos(theta+travel_offset) + pc(1),Rinner*sin(theta+travel_offset) + pc(2),'-','color',[0.7 0.7 0.7]);
plot(Rmin*cos(theta+travel_offset) + pc(1),Rmin*sin(theta+travel_offset) + pc(2),'r-');
plot(Rmax*cos(theta+travel_offset+theta_max_offset) + pc(1),Rmax*sin(theta+travel_offset+theta_max_offset) + pc(2),'b-');
plot(Router*cos(theta+travel_offset) + pc(1),Router*sin(theta+travel_offset) + pc(2),'-','color',[0.7 0.7 0.7]);
% plot(Rmax*cos(theta+travel_offset) + pc(1),Rmax*sin(theta+travel_offset) + pc(2),'b-');

% Now determine where the front corners of the vehicle are at each moment
plf = zeros(N,2);
prf = zeros(N,2);
plr = zeros(N,2);
prr = zeros(N,2);
plf(:,1) = pv(:,1) + vehicle.a*cos(theta+h0) + vehicle.d/2*cos(theta+h0+pi/2);
plf(:,2) = pv(:,2) + vehicle.a*sin(theta+h0) + vehicle.d/2*sin(theta+h0+pi/2);
prf(:,1) = pv(:,1) + vehicle.a*cos(theta+h0) + vehicle.d/2*cos(theta+h0-pi/2);
prf(:,2) = pv(:,2) + vehicle.a*sin(theta+h0) + vehicle.d/2*sin(theta+h0-pi/2);
plr(:,1) = pv(:,1) - vehicle.b*cos(theta+h0) + vehicle.d/2*cos(theta+h0+pi/2);
plr(:,2) = pv(:,2) - vehicle.b*sin(theta+h0) + vehicle.d/2*sin(theta+h0+pi/2);
prr(:,1) = pv(:,1) - vehicle.b*cos(theta+h0) + vehicle.d/2*cos(theta+h0-pi/2);
prr(:,2) = pv(:,2) - vehicle.b*sin(theta+h0) + vehicle.d/2*sin(theta+h0-pi/2);
% And plot the front corners continuously
% plot(plf(:,1),plf(:,2),'r-')
% plot(prf(:,1),prf(:,2),'b-')
% plot(plr(:,1),plr(:,2),'r-.')
% plot(prr(:,1),prr(:,2),'b-.')
% Plot a few of the rear corners
inds = 1;%;[1; floor(N/4); floor(N/2); floor(3*N/4)];%
for i = 1:length(inds)
    plot(pv(inds(i),1),pv(inds(i),2),'ko')
    plot([plf(inds(i),1) prf(inds(i),1) prr(inds(i),1) plr(inds(i),1) plf(inds(i),1)],...
        [plf(inds(i),2) prf(inds(i),2) prr(inds(i),2) plr(inds(i),2) plf(inds(i),2)],'k');
end

if 1 == getNewData
    pause();
    Nobstacles = 1;
    % Add an obstacle by clicking to generate the points and then turning the
    % points into a patch object
    obstacles = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});
    
    for obstacleInd = 1:Nobstacles
        [xobst,yobst] = ginput;
        obstacles(obstacleInd).pointsX = xobst;
        obstacles(obstacleInd).pointsY = yobst;
        obstacles(obstacleInd).color = [0.4 0.4 0.4];
    end
    % else
    %         xobst = [
    %             6.8899
    %             7.1416
    %             23.5009
    %             23.0695
    %             ];
    %         yobst = [
    %             0.9888
    %             -2.8943
    %             -3.8291
    %             8.9348
    %             ];
    
    %         xobst = [
    %             15.8899
    %             16.1416
    %             10000
    %             23.0695
    %             ];
    %         yobst = [
    %             0.9888
    %             -2.8943
    %             -3.8291
    %             8.9348
    %             ];
    %
    %         obstacles(1).pointsX = xobst;
    %         obstacles(1).pointsY = yobst;
    %         obstacles(1).color = [0.4 0.4 0.4];
end
% Plot the patch object over the trajectory lines using the patch plotting
% utility
[h,hpts] = fcn_Patch_plotPatch(obstacles,1);
axis auto

%% Check for a collision at the initial configuration
initial_collision_flag = 0;
carBoundBox = [plf(1,:); prf(1,:); prr(1,:); plr(1,:); plf(1,:)];

for obstacleInd = 1:Nobstacles
    [XI,YI] = polyxpoly(carBoundBox(:,1),carBoundBox(:,2),obstacles(obstacleInd).pointsX([1:end 1]),obstacles(obstacleInd).pointsY([1:end 1]));
    
    if ~isempty(XI)
        initial_collision_flag = 1;
        fprintf(1,'Vehicle initial position is in conflict with obstacle %d.\n',obstacleInd);
        collisionPoly = intersect(polyshape(carBoundBox(:,1),carBoundBox(:,2)),polyshape(obstacles(obstacleInd).pointsX([1:end 1]),obstacles(obstacleInd).pointsY([1:end 1])));
        plot(collisionPoly,'facecolor','red')
    end
end

%% Determine the nearest collision
if 0 == initial_collision_flag
    tic
    [collFlags,collTime,collAngle,collLoc,clearance,bodyCollLoc] = fcn_Patch_checkCollisions(x0,vehicle,obstacles);
    ET = toc;
    
    fprintf(1,'Determined collision geometry for %d objects in %0.3f seconds.\n',length(collFlags),ET);
    
    % Plot the car body
    Rabs = abs(R);
    forwardAngleOffset = sign(R)*pi/2;
    
    for collInd = 1:length(collFlags)
        bodyLoc(1,1) = Rabs*cos(collAngle(collInd)) + pc(1) + vehicle.a*cos(collAngle(collInd)+forwardAngleOffset) - vehicle.d/2*cos(collAngle(collInd));
        bodyLoc(1,2) = Rabs*sin(collAngle(collInd)) + pc(2) + vehicle.a*sin(collAngle(collInd)+forwardAngleOffset) - vehicle.d/2*sin(collAngle(collInd));
        bodyLoc(2,1) = Rabs*cos(collAngle(collInd)) + pc(1) + vehicle.a*cos(collAngle(collInd)+forwardAngleOffset) + vehicle.d/2*cos(collAngle(collInd));
        bodyLoc(2,2) = Rabs*sin(collAngle(collInd)) + pc(2) + vehicle.a*sin(collAngle(collInd)+forwardAngleOffset) + vehicle.d/2*sin(collAngle(collInd));
        bodyLoc(3,1) = Rabs*cos(collAngle(collInd)) + pc(1) - vehicle.b*cos(collAngle(collInd)+forwardAngleOffset) + vehicle.d/2*cos(collAngle(collInd));
        bodyLoc(3,2) = Rabs*sin(collAngle(collInd)) + pc(2) - vehicle.b*sin(collAngle(collInd)+forwardAngleOffset) + vehicle.d/2*sin(collAngle(collInd));
        bodyLoc(4,1) = Rabs*cos(collAngle(collInd)) + pc(1) - vehicle.b*cos(collAngle(collInd)+forwardAngleOffset) - vehicle.d/2*cos(collAngle(collInd));
        bodyLoc(4,2) = Rabs*sin(collAngle(collInd)) + pc(2) - vehicle.b*sin(collAngle(collInd)+forwardAngleOffset) - vehicle.d/2*sin(collAngle(collInd));
        plot(bodyLoc([1:end 1],1),bodyLoc([1:end 1],2),'b-','linewidth',1);
        
        % Plot the collision or closest clearance point
        plot(collLoc(collInd,1),collLoc(collInd,2),'r*')
        
        if 1 == collFlags(collInd)
            fprintf(1,'  For object %d, found a collision with TTC of %0.2f seconds at (%0.2f,%0.2f)\n',collInd,collTime(collInd),collLoc(collInd,1),collLoc(collInd,2));
        else
            fprintf(1,'  For object %d, no collision detected. Smallest clearance distance is %0.2f units at (%0.2f,%0.2f), occurring at %0.2f seconds.\n',collInd,clearance(collInd),collLoc(collInd,1),collLoc(collInd,2),collTime(collInd));
        end
    end
    
    [~,firstColl] = conditionalMin(collTime,collFlags,'== 1');
    
    % Create another copy of the figure
    figure(1)
    a1 = gca;
    f2 = figure(2);
    clf
    a2 = copyobj(a1,f2);
    figure(2)
    axis([collLoc(firstColl,1)-2 collLoc(firstColl,1)+2 collLoc(firstColl,2)-2 collLoc(firstColl,2) + 2]);
    % Create a legend
    %legend('Vehicle Start Point','Vehicle Trajectory','Inner Vehicle Bound','Outer Vehicle Bound','Vehicle CG','Vehicle Outline at Start','Obstacle','Obstacle Vertices','Vehicle Outline at Collision','Collision Point','location','best')
end