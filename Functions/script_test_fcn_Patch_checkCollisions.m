% Script to generate vehicle path and collision data using a constant
% radius circular model of the vehicle behavior. Currently assumes the
% vehicle CG travels on the circle and that there is zero sideslip angle.
% Could be extended to include vehicle sideslip, though the approach
% would have to use either max sideslip angle or change the bounds of the
% vehicle as it goes around.


% Remaining issues: clockwise travel

%clearvars

% Vehicle trajectory information
vx = 20;        % longitudinal speed (m/s)
R = 15;         % path radius (m) with sign (+ left, - right)
p0 = [0,0];     % initial position of vehicle (m,m)
h0 = 0;      % initial heading of vehicle (rad)
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
Rmax = sqrt((R+vehicle.d/2*sign(R))^2 + max(vehicle.a,vehicle.b)^2);
theta_max_offset = atan2(-vehicle.b,R+vehicle.d/2*sign(R));
pmin = zeros(N,2);
pmax = zeros(N,2);
pmin(:,1) = Rmin*cos(theta+h0-pi/2) + pc(1);
pmin(:,2) = Rmin*sin(theta+h0-pi/2) + pc(2);
pmax(:,1) = Rmax*cos(theta+h0-pi/2+theta_max_offset) + pc(1);
pmax(:,2) = Rmax*sin(theta+h0-pi/2+theta_max_offset) + pc(2);
plot(pmin(:,1),pmin(:,2),'r-');
plot(pmax(:,1),pmax(:,2),'b-');

% Now determine where the front corners of the vehicle are at each moment
plf = zeros(N,2);
prf = zeros(N,2);
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
inds = 1%;[1; floor(N/4); floor(N/2); floor(3*N/4)];
for i = 1:length(inds)
    plot(pv(inds(i),1),pv(inds(i),2),'ko')
    plot([plf(inds(i),1) prf(inds(i),1) prr(inds(i),1) plr(inds(i),1) plf(inds(i),1)],...
        [plf(inds(i),2) prf(inds(i),2) prr(inds(i),2) plr(inds(i),2) plf(inds(i),2)],'k');
end

% Add an obstacle by clicking to generate the points and then turning the
% points into a patch object
%axis([9 12 25 27])
%[xobst,yobst] = ginput;
obstacles = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});
obstacles(1).pointsX = xobst;
obstacles(1).pointsY = yobst;
obstacles(1).color = [0.4 0.4 0.4];
% Plot the patch object over the trajectory lines using the patch plotting
% utility
[h,hpts] = fcn_Patch_plotPatch(obstacles,1);

%% Determine the nearest collision
[collTime,collAngle,collLoc,~] = fcn_Patch_checkCollisions(x0,vehicle,obstacles);

% Plot the car body
bodyLoc(1,1) = R*cos(collAngle) + pc(1) + vehicle.a*cos(collAngle+pi/2) - vehicle.d/2*cos(collAngle);
bodyLoc(1,2) = R*sin(collAngle) + pc(2) + vehicle.a*sin(collAngle+pi/2) - vehicle.d/2*sin(collAngle);
bodyLoc(2,1) = R*cos(collAngle) + pc(1) + vehicle.a*cos(collAngle+pi/2) + vehicle.d/2*cos(collAngle);
bodyLoc(2,2) = R*sin(collAngle) + pc(2) + vehicle.a*sin(collAngle+pi/2) + vehicle.d/2*sin(collAngle);
bodyLoc(3,1) = R*cos(collAngle) + pc(1) - vehicle.b*cos(collAngle+pi/2) + vehicle.d/2*cos(collAngle);
bodyLoc(3,2) = R*sin(collAngle) + pc(2) - vehicle.b*sin(collAngle+pi/2) + vehicle.d/2*sin(collAngle);
bodyLoc(4,1) = R*cos(collAngle) + pc(1) - vehicle.b*cos(collAngle+pi/2) - vehicle.d/2*cos(collAngle);
bodyLoc(4,2) = R*sin(collAngle) + pc(2) - vehicle.b*sin(collAngle+pi/2) - vehicle.d/2*sin(collAngle);
plot(bodyLoc([1:end 1],1),bodyLoc([1:end 1],2),'b-','linewidth',1);

% Plot the collision point
plot(collLoc(1),collLoc(2),'r*')

% Create a legend
legend('Vehicle Start Point','Vehicle Trajectory','Inner Vehicle Bound','Outer Vehicle Bound','Vehicle CG','Vehicle Outline at Start','Obstacle','Obstacle Vertices','Vehicle Outline at Collision','Collision Point','location','best')
