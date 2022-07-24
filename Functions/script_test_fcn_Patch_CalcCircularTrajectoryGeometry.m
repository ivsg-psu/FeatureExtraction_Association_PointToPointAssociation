% Script to exercise the function fcn_Patch_CalcCircularTrajectoryGeometry

clearvars

% Vehicle trajectory information
vx = 20;        % longitudinal speed (m/s)
R = 20;         % path radius (m) with sign (+ left, - right)
p0 = [10,-5];     % initial position of vehicle (m,m)
h0 = pi/2;     % initial heading of vehicle (rad)
a0 = 5*pi/180; % vehicle body slip angle (rad)
tf = 1;         % time horizon to check (s)
% Create a vector of the trajectory information
x0 = [p0'; h0; a0; vx; R];
% Vehicle dimensional information
vehicle.dr = 2.2;       % CG-front bumper distance (m)
vehicle.df = 3.5;        % CG-rear bumper distance (m)
vehicle.w = 2.3;        % vehicle width (m)

% Set up a new figure (or clear the existing one)
figure(1)
clf
hold on
grid on
axis equal

% Plot the vehicle CG at the initial location
plot(p0(1),p0(2),'ko','markersize',8)
% Plot a plus sign over the circle to make something akin to a CG symbol
plot(p0(1),p0(2),'k+','markersize',8)

% Calculate center point of vehicle trajectory circle
pc(1) = p0(1) + R*cos(h0+pi/2);
pc(2) = p0(2) + R*sin(h0+pi/2);

% Generate a trajectory for the given time horizon
del_t = 0.01;
t = (0:del_t:tf)';
theta = vx*t/R + h0-pi/2;
N = length(t);
pv = zeros(N,2);
pv(:,1) = R*cos(theta) + pc(1);
pv(:,2) = R*sin(theta) + pc(2);
% Plot the trajectory
plot(pv(:,1),pv(:,2),'k-.')

% Calculate the various pertinent radii and corner points of the vehicle
[radii,vehicleBB,radiiFlags] = fcn_Patch_CalcCircularTrajectoryGeometry(x0,vehicle);

% Parse out which radii are which
Rinside = sign(R)*radii(1);
RoutsideFront = sign(R)*radii(2);
RoutsideRear = sign(R)*radii(3);
Rmin = sign(R)*radii(6);
Rmax = sign(R)*radii(7);

% Calculate the offset to the angular position based on the back end of the
% vehicle (in order to plot the extra portion of the clearance curves)
rear_offset = 1.1*sign(R)*atan2(-vehicle.dr,sign(R)*R+vehicle.w/2);
% Create an angular range over which to plot the various radii
theta_arcs = linspace(rear_offset,theta(end),100);
% Plot the various radii 
plot(Rinside*cos(theta_arcs) + pc(1),Rinside*sin(theta_arcs) + pc(2),'-','color',[0.7 0.7 0.7]);
plot(RoutsideFront*cos(theta_arcs) + pc(1),RoutsideFront*sin(theta_arcs) + pc(2),'-','color',[0.7 0.7 0.7]);
plot(RoutsideRear*cos(theta_arcs) + pc(1),RoutsideRear*sin(theta_arcs) + pc(2),'-','color',[0.7 0.7 0.7]);
plot(Rmin*cos(theta_arcs) + pc(1),Rmin*sin(theta_arcs) + pc(2),'r-');
plot(Rmax*cos(theta_arcs) + pc(1),Rmax*sin(theta_arcs) + pc(2),'b-');

% Finally, plot the vehicle body over the top of it all
plot(vehicleBB([1:4 1],1),vehicleBB([1:4 1],2),'k.-')
% If the inner tangent is the limiting inner radii, plot the tangent point
if 5 == radiiFlags(1)
    plot(vehicleBB(5,1),vehicleBB(5,2),'r*')
end