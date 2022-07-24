function fcn_Patch_plotCircularCollisions(x0,vehicle,patchArray,collAngle,collLoc,t_f,varargin)
% fcn_Patch_plotPatch
% Plots a visual representation of the objects in a patch structure array
%
% FORMAT: 
%
%       fcn_Patch_plotStraightCollisions(,(fig_num))
%
% INPUTS:
%
%      collFlags: an N x 1 vector of flags denoting whether there is a
%           collision with each of the objects
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
%       See the script: script_test_fcn_Patch_plotCircularCollisions.m for a full test
%       suite. 
%
% This function was written by Shashank Garikipati

% Revision history:
%     2022_07_12
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
    if nargin < 6 || nargin > 7
        error('Incorrect number of input arguments')
    end
    
end

% Did the user provide a figure number?
if 7 == nargin
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

p0 = x0(1:2)
CoG0 = x0(3);
a0 = x0(4);
vx = x0(5);
R = x0(6);

pc(1) = p0(1) + R*cos(CoG0+pi/2);
pc(2) = p0(2) + R*sin(CoG0+pi/2);

del_t = 0.01;
t = (0:del_t:t_f)';
theta = vx*t/R + CoG0 - sign(R)*pi/2;

% Iterate over all of objects
for obstacleInd = 1:size(patchArray,2)
    % Determine the vehicle heading at the point of collision or near miss
    % using the angle of the initial vehicle position, theta(1)
    collHeading = theta(1) + sign(R)*collAngle(obstacleInd) + (pi/2+a0);
    % Determine the CG location at the point of collision or near miss
    collCG = [abs(R)*cos(theta(1) + sign(R)*collAngle(obstacleInd)) ...
              abs(R)*sin(theta(1) + sign(R)*collAngle(obstacleInd))] + pc;
    % Plot the vehicle body in the position of the collision or near miss
    plotVehicleBB(collCG,collHeading,vehicle,R,fig_num)
    
    % Plot the collision or closest clearance point
    plot(collLoc(obstacleInd,1),collLoc(obstacleInd,2),'r*')
end

end

function plotVehicleBB(pCG,heading,vehicle,R,figHandle)

% Enumerate the body positions and directions
LF = 1; RF = 2; RR = 3; LR = 4; X = 1; Y = 2;

% Address the specified figure
figure(figHandle)
% Plot the CG location
plot(pCG(X),pCG(Y),'ko')
plot(pCG(X),pCG(Y),'k+')

if R < 0
    heading = heading + pi;
end

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
% Plot the vehicle bounding box
plot(vehicleBB([1:end 1],X),vehicleBB([1:end 1],Y),'b-','linewidth',1);
end



