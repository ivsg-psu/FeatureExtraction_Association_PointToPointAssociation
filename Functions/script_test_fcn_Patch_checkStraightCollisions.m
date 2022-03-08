% Script to generate vehicle path and collision data using a straight line
% model of the vehicle behavior. Currently assumes that there is zero
% sideslip angle.

% Notation: px denotes points. For example, p0 is initial point, pc is the
% center point, etc.

getNewData = 1;
if 1 == getNewData
    clearvars
    getNewData = 1;
end

% Vehicle trajectory information
vx = 20;        % longitudinal speed (m/s)
p0 = [2,8];     % initial position of vehicle (m,m)
h0 = 3*pi/5;      % initial heading of vehicle (rad)
tf = 4;         % time horizon to check (s)
% Create a vector of the trajectory information
x0 = [p0'; h0; vx];
% Vehicle dimensional information
vehicle.df = 1.8;       % CG-front bumper distance (m)
vehicle.dr = 2.2;        % CG-rear bumper distance (m)
vehicle.w = 2.0;        % vehicle width (m)

% Plot the points
figure(1)
clf
hold on
grid on
axis equal
plot(p0(1),p0(2),'k*')

% Generate a trajectory for the given time horizon
del_t = 0.01;
t = (0:del_t:tf)';
d = vx*t;
N = length(t);
pv = zeros(N,2);
pv(:,1) = d*cos(h0) + p0(1);
pv(:,2) = d*sin(h0) + p0(2);
% And plot the trajectory
plot(pv(:,1),pv(:,2),'k-.')

% Plot the outer limits of the vehicle trajectory
%plot(pv(:,1)+vehicle.w/2*sin(h0),pv(:,2)-vehicle.w/2*cos(h0),'b');
plot(pv(:,1)-vehicle.dr*cos(h0)-vehicle.w/2*sin(h0),pv(:,2)-vehicle.dr*sin(h0)+vehicle.w/2*cos(h0),'b');
plot(pv(:,1)-vehicle.dr*cos(h0)+vehicle.w/2*sin(h0),pv(:,2)-vehicle.dr*sin(h0)-vehicle.w/2*cos(h0),'r');
%axis([-10+p0(1) 10+p0(1) -10+p0(2) 10+p0(2)])

% Now determine where the front corners of the vehicle are at each moment
% plf = zeros(N,2);
% prf = zeros(N,2);
% plr = zeros(N,2);
% prr = zeros(N,2);
% plf(:,1) = pv(:,1) + vehicle.df*cos(theta+h0) + vehicle.w/2*cos(theta+h0+pi/2);
% plf(:,2) = pv(:,2) + vehicle.df*sin(theta+h0) + vehicle.w/2*sin(theta+h0+pi/2);
% prf(:,1) = pv(:,1) + vehicle.df*cos(theta+h0) + vehicle.w/2*cos(theta+h0-pi/2);
% prf(:,2) = pv(:,2) + vehicle.df*sin(theta+h0) + vehicle.w/2*sin(theta+h0-pi/2);
% plr(:,1) = pv(:,1) - vehicle.dr*cos(theta+h0) + vehicle.w/2*cos(theta+h0+pi/2);
% plr(:,2) = pv(:,2) - vehicle.dr*sin(theta+h0) + vehicle.w/2*sin(theta+h0+pi/2);
% prr(:,1) = pv(:,1) - vehicle.dr*cos(theta+h0) + vehicle.w/2*cos(theta+h0-pi/2);
% prr(:,2) = pv(:,2) - vehicle.dr*sin(theta+h0) + vehicle.w/2*sin(theta+h0-pi/2);
% And plot the front corners continuously
% plot(plf(:,1),plf(:,2),'r-')
% plot(prf(:,1),prf(:,2),'b-')
% plot(plr(:,1),plr(:,2),'r-.')
% plot(prr(:,1),prr(:,2),'b-.')
% Plot a few of the rear corners
% inds = 1;%;[1; floor(N/4); floor(N/2); floor(3*N/4)];%
% for i = 1:length(inds)
%     plot(pv(inds(i),1),pv(inds(i),2),'ko')
%     plot([plf(inds(i),1) prf(inds(i),1) prr(inds(i),1) plr(inds(i),1) plf(inds(i),1)],...
%         [plf(inds(i),2) prf(inds(i),2) prr(inds(i),2) plr(inds(i),2) plf(inds(i),2)],'k');
% end

if 1 == getNewData
    Nobstacles = 2;
    % Add an obstacle by clicking to generate the points and then turning the
    % points into a patch object
    obstacles = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});
    
    for obstacleInd = 1:Nobstacles
        [xobst,yobst] = ginput;
        obstacles(obstacleInd).pointsX = xobst;
        obstacles(obstacleInd).pointsY = yobst;
        obstacles(obstacleInd).color = [0.4 0.4 0.4];
    end
end
% Plot the patch object over the trajectory lines using the patch plotting
% utility
[h,hpts] = fcn_Patch_plotPatch(obstacles,1);
axis auto

%% Check for a collision at the initial configuration
initial_collision_flag = 0;
% matA = [vehicle.df*ones(2,1) vehicle.w/2*[1; -1]; -vehicle.dr*ones(2,1) vehicle.w/2*[-1; 1]; vehicle.df vehicle.w/2];
% matB = [cosd(h0) sind(h0); sind(h0) -cosd(h0)];
% carBoundBox = matA*matB;
% 
% for obstacleInd = 1:Nobstacles
%     [XI,YI] = polyxpoly(carBoundBox(:,1),carBoundBox(:,2),obstacles(obstacleInd).pointsX([1:end 1]),obstacles(obstacleInd).pointsY([1:end 1]));
%     
%     if ~isempty(XI)
%         initial_collision_flag = 1;
%         fprintf(1,'Vehicle initial position is in conflict with obstacle %d.\n',obstacleInd);
%         collisionPoly = intersect(polyshape(carBoundBox(:,1),carBoundBox(:,2)),polyshape(obstacles(obstacleInd).pointsX([1:end 1]),obstacles(obstacleInd).pointsY([1:end 1])));
%         plot(collisionPoly,'facecolor','red')
%     end
% end

%% Determine the nearest collision
if 0 == initial_collision_flag
    tic
    [collFlags,collTime,collLoc,clearance,bodyCollLoc] = fcn_Patch_checkStraightCollisions(x0,vehicle,obstacles);
    ET = toc;
    
    fprintf(1,'Determined collision geometry for %d objects in %0.3f seconds.\n',length(collFlags),ET);
    
%plot(pv(:,1)-vehicle.dr*cos(h0)-vehicle.w/2*sin(h0),pv(:,2)-vehicle.dr*sin(h0)+vehicle.w/2*cos(h0),'b');
%plot(pv(:,1)-vehicle.dr*cos(h0)+vehicle.w/2*sin(h0),pv(:,2)-vehicle.dr*sin(h0)-vehicle.w/2*cos(h0),'r');
    % Determine the offsets from the collision (or near miss) location in
    % order to plot the vehicle body in the correct location
    % Plot the car body
    for collInd = 1:length(collFlags)
        plotOffset(1) = bodyCollLoc(collInd,1);
        if isnan(clearance(collInd))
            plotOffset(2) = bodyCollLoc(collInd,2);
        else
            plotOffset(2) = bodyCollLoc(collInd,2) + sign(bodyCollLoc(collInd,2))*clearance(collInd);
        end
        bodyLoc(1,1) = collLoc(collInd,1) + (vehicle.df - plotOffset(1))*cos(h0) - (vehicle.w/2-plotOffset(2))*sin(h0);
        bodyLoc(1,2) = collLoc(collInd,2) + (vehicle.df - plotOffset(1))*sin(h0) + (vehicle.w/2-plotOffset(2))*cos(h0);
        bodyLoc(2,1) = collLoc(collInd,1) + (vehicle.df - plotOffset(1))*cos(h0) - (-vehicle.w/2-plotOffset(2))*sin(h0);
        bodyLoc(2,2) = collLoc(collInd,2) + (vehicle.df - plotOffset(1))*sin(h0) + (-vehicle.w/2-plotOffset(2))*cos(h0);
        bodyLoc(3,1) = collLoc(collInd,1) + (-vehicle.dr - plotOffset(1))*cos(h0) - (-vehicle.w/2-plotOffset(2))*sin(h0);
        bodyLoc(3,2) = collLoc(collInd,2) + (-vehicle.dr - plotOffset(1))*sin(h0) + (-vehicle.w/2-plotOffset(2))*cos(h0);
        bodyLoc(4,1) = collLoc(collInd,1) + (-vehicle.dr - plotOffset(1))*cos(h0) - (vehicle.w/2-plotOffset(2))*sin(h0);
        bodyLoc(4,2) = collLoc(collInd,2) + (-vehicle.dr - plotOffset(1))*sin(h0) + (vehicle.w/2-plotOffset(2))*cos(h0);
        figure(1)
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
%     % Create a legend
%     %legend('Vehicle Start Point','Vehicle Trajectory','Inner Vehicle Bound','Outer Vehicle Bound','Vehicle CG','Vehicle Outline at Start','Obstacle','Obstacle Vertices','Vehicle Outline at Collision','Collision Point','location','best')
end