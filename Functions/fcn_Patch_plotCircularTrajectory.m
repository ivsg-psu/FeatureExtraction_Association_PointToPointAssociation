xfunction fcn_Patch_plotCircularTrajectory(t_f,x0,vehicle,varargin)
% fcn_Patch_plotPatch
% Plots a visual representation of the objects in a patch structure array
%
% FORMAT: 
%
%       fcn_Patch_plotCircularTrajectory(t_h,x0,vehicle,(fig_num)
%
% INPUTS:
%
%      tf: [1x1] scalar, representing time horizon to check
%      x0: [4x1] vector of [p0',h0,vx]
%         vx: [1x1] scalar, representing longitudinal speed
%         p0: [1x1] scalar, representing initial position of the vehicle
%         h0: [1x1] scalar, representing initial heading of the vehicle
%      vehicle: [1x1] cell struct with vehicle properties
%
%      (OPTIONAL INPUTS)
%      fig_num: a figure number to plot into
%
% OUTPUTS:
%
%
% DEPENDENCIES:
%
%
% EXAMPLES:
%      
%       See the script: script_test_fcn_Patch_plotStraightTragectory.m for a full test
%       suite. 
%
% This function was written by Shashank Garikipati

% Revision history:
%     2022_07_06
%     -- wrote the code

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
    if nargin < 3 || nargin > 4
        error('Incorrect number of input arguments')
    end
    
end

% Did the user provide a figure number?
if 3 < nargin
    fig_num = varargin{1};
    figure(fig_num);
else    
    fig = figure;
    fig_num = fig.Number;
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
 
figure(fig_num);
grid on
hold on
axis equal

% Enumerate the body positions and directions
LF = 1; RF = 2; RR = 3; LR = 4; X = 1; Y = 2;

p0 = x0(1:2)';
CoG0 = x0(3);
a0 = x0(4);
vx = x0(5);
R = x0(6);

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

% Plot a plus sign over the circle to make something akin to a CG symbol
plot(p0(1),p0(2),'k+','markersize',8)

% Calculate center point of vehicle trajectory circle
pc(1) = p0(1) + R*cos(CoG0+pi/2);
pc(2) = p0(2) + R*sin(CoG0+pi/2);

% Generate a circular trajectory for the given time horizon
del_t = 0.01;
t = (0:del_t:t_f)';
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

plotVehicleBB(p0,CoG0+a0,vehicle,fig_num)
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

body1Loc(1,1) = pCG(X) + (vehicle.df - 0)*cos(heading) - (((vehicle.w/2-(1/2)*vehicle.w/2))-0)*sin(heading);
body1Loc(1,2) = pCG(Y) + (vehicle.df - 0)*sin(heading) + (((vehicle.w/2-(1/2)*vehicle.w/2))-0)*cos(heading);
body1Loc(2,1) = pCG(X) + ((vehicle.df-(1/3.5)*vehicle.df) - 0)*cos(heading) - (vehicle.w/2-0)*sin(heading);
body1Loc(2,2) = pCG(Y) + ((vehicle.df-(1/3.5)*vehicle.df) - 0)*sin(heading) + (vehicle.w/2-0)*cos(heading);
body1Loc(3,1) = pCG(X) + ((-vehicle.dr/2-(1/8)*vehicle.dr/2) - 0)*cos(heading) - (vehicle.w/2-0)*sin(heading);
body1Loc(3,2) = pCG(Y) + ((-vehicle.dr/2-(1/8)*vehicle.dr/2) - 0)*sin(heading) + (vehicle.w/2-0)*cos(heading);
body1Loc(4,1) = pCG(X) + (-vehicle.dr - 0)*cos(heading) - ((vehicle.w/2-(1/2)*vehicle.w/2)-0)*sin(heading);
body1Loc(4,2) = pCG(Y) + (-vehicle.dr - 0)*sin(heading) + ((vehicle.w/2-(1/2)*vehicle.w/2)-0)*cos(heading);

body1Loc(5,1) = pCG(X) + (-vehicle.dr - 0)*cos(heading) - (((-vehicle.w/2+(1/2)*vehicle.w/2))-0)*sin(heading);
body1Loc(5,2) = pCG(Y) + (-vehicle.dr - 0)*sin(heading) + (((-vehicle.w/2+(1/2)*vehicle.w/2))-0)*cos(heading);
body1Loc(6,1) = pCG(X) + ((-vehicle.dr/2-(1/8)*vehicle.dr/2) - 0)*cos(heading) - (-vehicle.w/2-0)*sin(heading);
body1Loc(6,2) = pCG(Y) + ((-vehicle.dr/2-(1/8)*vehicle.dr/2) - 0)*sin(heading) + (-vehicle.w/2-0)*cos(heading);
body1Loc(7,1) = pCG(X) + ((vehicle.df-(1/3.5)*vehicle.df) - 0)*cos(heading) - (-vehicle.w/2-0)*sin(heading);
body1Loc(7,2) = pCG(Y) + ((vehicle.df-(1/3.5)*vehicle.df) - 0)*sin(heading) + (-vehicle.w/2-0)*cos(heading);
body1Loc(8,1) = pCG(X) + (vehicle.df - 0)*cos(heading) - (((-vehicle.w/2+(1/2)*vehicle.w/2))-0)*sin(heading);
body1Loc(8,2) = pCG(Y) + (vehicle.df - 0)*sin(heading) + (((-vehicle.w/2+(1/2)*vehicle.w/2))-0)*cos(heading);
plot(body1Loc([1:end 1],1),body1Loc([1:end 1],2),'k','linewidth',1);


body2Loc(1,1) = pCG(X) + ((vehicle.df-(1/3)*vehicle.df) - 0)*cos(heading) - (vehicle.w/2-0)*sin(heading);
body2Loc(1,2) = pCG(Y) + ((vehicle.df-(1/3)*vehicle.df) - 0)*sin(heading) + (vehicle.w/2-0)*cos(heading);
body2Loc(2,1) = pCG(X) + ((vehicle.df-(1/3.8)*vehicle.df) - 0)*cos(heading) - ((vehicle.w/2-(1/2)*vehicle.w/2)-0)*sin(heading);
body2Loc(2,2) = pCG(Y) + ((vehicle.df-(1/3.8)*vehicle.df) - 0)*sin(heading) + ((vehicle.w/2-(1/2)*vehicle.w/2)-0)*cos(heading);
body2Loc(3,1) = pCG(X) + ((vehicle.df-(1/3.8)*vehicle.df) - 0)*cos(heading) - ((-vehicle.w/2+(1/2)*vehicle.w/2)-0)*sin(heading);
body2Loc(3,2) = pCG(Y) + ((vehicle.df-(1/3.8)*vehicle.df) - 0)*sin(heading) + ((-vehicle.w/2+(1/2)*vehicle.w/2)-0)*cos(heading);
body2Loc(4,1) = pCG(X) + ((vehicle.df-(1/3)*vehicle.df) - 0)*cos(heading) - ((-vehicle.w/2)-0)*sin(heading);
body2Loc(4,2) = pCG(Y) + ((vehicle.df-(1/3)*vehicle.df) - 0)*sin(heading) + ((-vehicle.w/2)-0)*cos(heading);

body2Loc(5,1) = pCG(X) + ((vehicle.df-(1/1.8)*vehicle.df) - 0)*cos(heading) - ((-vehicle.w/2+(1/4)*vehicle.w/2)-0)*sin(heading);
body2Loc(5,2) = pCG(Y) + ((vehicle.df-(1/1.8)*vehicle.df) - 0)*sin(heading) + ((-vehicle.w/2+(1/4)*vehicle.w/2)-0)*cos(heading);
body2Loc(6,1) = pCG(X) + ((vehicle.df-(1/1.8)*vehicle.df) - 0)*cos(heading) - ((vehicle.w/2-(1/4)*vehicle.w/2)-0)*sin(heading);
body2Loc(6,2) = pCG(Y) + ((vehicle.df-(1/1.8)*vehicle.df) - 0)*sin(heading) + ((vehicle.w/2-(1/4)*vehicle.w/2)-0)*cos(heading);

plot(body2Loc([1:end 1],1),body2Loc([1:end 1],2),'k','linewidth',1);

body3Loc(1,1) = pCG(X) + ((-vehicle.dr/6) - 0)*cos(heading) - ((vehicle.w/2-(1.5/4)*vehicle.w/2)-0)*sin(heading);
body3Loc(1,2) = pCG(Y) + ((-vehicle.dr/6) - 0)*sin(heading) + ((vehicle.w/2-(1.5/4)*vehicle.w/2)-0)*cos(heading);
body3Loc(2,1) = pCG(X) + ((-vehicle.dr/6) - 0)*cos(heading) - (-(vehicle.w/2-(1.5/4)*vehicle.w/2)-0)*sin(heading);
body3Loc(2,2) = pCG(Y) + ((-vehicle.dr/6) - 0)*sin(heading) + (-(vehicle.w/2-(1.5/4)*vehicle.w/2)-0)*cos(heading);
body3Loc(3,1) = pCG(X) + ((-vehicle.dr/1.5) - 0)*cos(heading) - (-(vehicle.w/2-(1.5/4)*vehicle.w/2)-0)*sin(heading);
body3Loc(3,2) = pCG(Y) + ((-vehicle.dr/1.5) - 0)*sin(heading) + (-(vehicle.w/2-(1.5/4)*vehicle.w/2)-0)*cos(heading);
body3Loc(4,1) = pCG(X) + ((-vehicle.dr/1.4) - 0)*cos(heading) - (-0)*sin(heading);
body3Loc(4,2) = pCG(Y) + ((-vehicle.dr/1.4) - 0)*sin(heading) + (-0)*cos(heading);
body3Loc(5,1) = pCG(X) + ((-vehicle.dr/1.5) - 0)*cos(heading) - ((vehicle.w/2-(1.5/4)*vehicle.w/2)-0)*sin(heading);
body3Loc(5,2) = pCG(Y) + ((-vehicle.dr/1.5) - 0)*sin(heading) + ((vehicle.w/2-(1.5/4)*vehicle.w/2)-0)*cos(heading);

plot(body3Loc([1:end 1],1),body3Loc([1:end 1],2),'k','linewidth',1);
end





