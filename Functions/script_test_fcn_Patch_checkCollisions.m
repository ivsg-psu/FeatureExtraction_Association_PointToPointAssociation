% Script to generate vehicle path and collision data using a constant
% radius circular model of the vehicle behavior. Currently assumes the
% vehicle CG travels on the circle and that there is zero sideslip angle.
% Could be extended to include vehicle sideslip, though the approach
% would have to use either max sideslip angle or change the bounds of the
% vehicle as it goes around.


% Remaining issues: misses edges that cross the outside of the trajectory
% but have no vertices inside (see screen shot) because the line-circle
% intersection function rejects edges with two intersections. The routine
% also does not (yet) find the minimum clearance for misses

% Update: seems to incorrectly get the closest point on an edge for a
% near miss

clearvars

% Vehicle trajectory information
vx = 20;        % longitudinal speed (m/s)
R = 15;         % path radius (m) with sign (+ left, - right)
p0 = [0,0];     % initial position of vehicle (m,m)
h0 = -pi/4;      % initial heading of vehicle (rad)
tf = 3;         % time horizon to check (s)
% Create a vector of the trajectory information
x0 = [p0'; h0; vx; R];
% Vehicle dimensional information
vehicle.a = 1.75;       % CG-front axle distance (m)
vehicle.b = 2.5;        % CG-rear axle distance (m)
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
Rmax = sign(R)*sqrt((R+vehicle.d/2*sign(R))^2 + max(vehicle.a,vehicle.b)^2);
Rinner = sign(R)*sqrt((R-vehicle.d/2*sign(R))^2 + vehicle.a^2);
Router = sign(R)*sqrt((R+vehicle.d/2*sign(R))^2 + vehicle.a^2);
travel_offset = h0 - pi/2;
theta_max_offset = sign(R)*atan2(-vehicle.b,sign(R)*R+vehicle.d/2);
plot(Rinner*cos(theta+travel_offset) + pc(1),Rinner*sin(theta+travel_offset) + pc(2),'-','color',[0.7 0.7 0.7]);
plot(Rmin*cos(theta+travel_offset) + pc(1),Rmin*sin(theta+travel_offset) + pc(2),'r-');
plot(Rmax*cos(theta+travel_offset+theta_max_offset) + pc(1),Rmax*sin(theta+travel_offset+theta_max_offset) + pc(2),'b-');
plot(Router*cos(theta+travel_offset) + pc(1),Router*sin(theta+travel_offset) + pc(2),'-','color',[0.7 0.7 0.7]);

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

% Add an obstacle by clicking to generate the points and then turning the
% points into a patch object
[xobst,yobst] = ginput;
obstacles = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});
obstacles(1).pointsX = xobst;
obstacles(1).pointsY = yobst;
obstacles(1).color = [0.4 0.4 0.4];
% Plot the patch object over the trajectory lines using the patch plotting
% utility
[h,hpts] = fcn_Patch_plotPatch(obstacles,1);
axis auto

%% Determine the nearest collision
[collFlags,collTime,collAngle,collLoc,clearance] = fcn_Patch_checkCollisions(x0,vehicle,obstacles);

% Plot the car body
Rabs = abs(R);
forwardAngleOffset = sign(R)*pi/2;
bodyLoc(1,1) = Rabs*cos(collAngle) + pc(1) + vehicle.a*cos(collAngle+forwardAngleOffset) - vehicle.d/2*cos(collAngle);
bodyLoc(1,2) = Rabs*sin(collAngle) + pc(2) + vehicle.a*sin(collAngle+forwardAngleOffset) - vehicle.d/2*sin(collAngle);
bodyLoc(2,1) = Rabs*cos(collAngle) + pc(1) + vehicle.a*cos(collAngle+forwardAngleOffset) + vehicle.d/2*cos(collAngle);
bodyLoc(2,2) = Rabs*sin(collAngle) + pc(2) + vehicle.a*sin(collAngle+forwardAngleOffset) + vehicle.d/2*sin(collAngle);
bodyLoc(3,1) = Rabs*cos(collAngle) + pc(1) - vehicle.b*cos(collAngle+forwardAngleOffset) + vehicle.d/2*cos(collAngle);
bodyLoc(3,2) = Rabs*sin(collAngle) + pc(2) - vehicle.b*sin(collAngle+forwardAngleOffset) + vehicle.d/2*sin(collAngle);
bodyLoc(4,1) = Rabs*cos(collAngle) + pc(1) - vehicle.b*cos(collAngle+forwardAngleOffset) - vehicle.d/2*cos(collAngle);
bodyLoc(4,2) = Rabs*sin(collAngle) + pc(2) - vehicle.b*sin(collAngle+forwardAngleOffset) - vehicle.d/2*sin(collAngle);
plot(bodyLoc([1:end 1],1),bodyLoc([1:end 1],2),'b-','linewidth',1);

% Plot the collision or closest clearance point
plot(collLoc(1),collLoc(2),'r*')

% Create another copy of the figure
figure(1)
a1 = gca;
f2 = figure(2);
clf
a2 = copyobj(a1,f2);
figure(2)
axis([collLoc(1)-2 collLoc(1)+2 collLoc(2)-2 collLoc(2) + 2]);
% Create a legend
%legend('Vehicle Start Point','Vehicle Trajectory','Inner Vehicle Bound','Outer Vehicle Bound','Vehicle CG','Vehicle Outline at Start','Obstacle','Obstacle Vertices','Vehicle Outline at Collision','Collision Point','location','best')

if 1 == collFlags(1)
    fprintf('Found a collision with TTC of %0.2f seconds at (%0.2f,%0.2f)\n',collTime,collLoc(1),collLoc(2));
else
    fprintf('No collision detected. Smallest clearance distance is %0.2f units at (%0.2f,%0.2f), occurring at %0.2f seconds.\n',clearance,collLoc(1),collLoc(2),collTime);
end