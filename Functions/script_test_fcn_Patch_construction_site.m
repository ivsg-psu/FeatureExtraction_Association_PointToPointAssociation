%%
%script_test_fcn_Patch_construction_site

%% Example 1
% Construction Zone with channeling device

fig_num = 1;
hold on

blockedLaneMarker = [0 11.1; 367 11.1];       % Lane line to represent blocked lane
plot(blockedLaneMarker(:,1),blockedLaneMarker(:,2),'k')
axis([0 367 -150 150])

topLaneMarker = [0 7.5; 367 7.5];           % Top lane marker
plot(topLaneMarker(:,1),topLaneMarker(:,2),'k--')

middleLaneMarker = [0 3.8; 367 3.8];        % Middle lane marker
plot(middleLaneMarker(:,1),middleLaneMarker(:,2),'k--')

bottomLaneMarker = [0 0.1; 367 0.1]         % Bottom lane marker
plot(bottomLaneMarker(:,1),bottomLaneMarker(:,2),'k')


% Vehicle trajectory information
vx = 25;        % longitudinal speed (m/s)
p0 = [0,5.55];     % initial position of vehicle (m,m)
h0 = 0;      % initial heading of vehicle (rad)
t_h = 14.7;         % time horizon to check (s)
% Create a vector of the trajectory information
x0 = [p0'; h0; vx];

% Vehicle dimensional information
vehicle.dr = 1.9;       % CG-front bumper distance (m)
vehicle.df = 3.0;        % CG-rear bumper distance (m)
vehicle.w = 1.8;        % vehicle width (m)

% First points for row of channeling device
barrels(1,1) = 50;barrels(1,2) = 7.5;
for i = 1:17    
    barrels(i+1,1) = barrels(i) + 18.5;
    barrels(i+1,2) = 7.5;
end

% Second points for  row of channeling device
barrels(19,1) = 0; barrels(19,2) = 11.1;
for i = 19:35
    barrels(i+1,1) = barrels(i) + 18.5;
    barrels(i+1,2) = 11.1;
end

% left diagonal points of channeling device
for i = 36:40
    barrels(i+1,1) = barrels(19,1) + (i-35)*8.3;
    barrels(i+1,2) = barrels(19,2) - (i-35)*0.55;
end

% right diagonal points of channeling device
for i = 41:45
    barrels(i+1,1) = barrels(36,1) + (i-40)*8.3;
    barrels(i+1,2) = barrels(36,2) - (i-40)*0.55;
end

% randomized orientation of the barrels
for j = 1:46
    angles(j,1) = 2*pi*rand
end

% converting points to barrels for channeling device
barrelsArray = fcn_Patch_fillBarrelAtPoint(barrels,angles);

% detect collision and show visualization
[collFlags,collTime,collLoc,clearance,bodyCollLoc] = fcn_Patch_checkStraightCollisions(x0,vehicle,barrelsArray,t_h,fig_num);

%% Example 2
%

fig_num2 = 2;
figure(fig_num2)
hold on

blockedLaneMarker = [0 11.1; 354 11.1];       % Lane line to represent blocked lane
plot(blockedLaneMarker(:,1),blockedLaneMarker(:,2),'k')
axis([0 367 -150 150])

topLaneMarker = [0 7.5; 354 7.5];           % Top lane marker
plot(topLaneMarker(:,1),topLaneMarker(:,2),'k--')

middleLaneMarker = [0 3.8; 354 3.8];        % Middle lane marker
plot(middleLaneMarker(:,1),middleLaneMarker(:,2),'k--')

bottomLaneMarker = [0 0.1; 354 0.1]         % Bottom lane marker
plot(bottomLaneMarker(:,1),bottomLaneMarker(:,2),'k')

% First points for row of channeling device
barrels2(1,1) = 45.72;barrels2(1,2) = 7.5;
for i = 1:17    
    barrels2(i+1,1) = barrels2(i) + 18.5;
    barrels2(i+1,2) = 7.5;
end


% randomized orientation of the barrels
for j = 1:46
    angles(j,1) = 2*pi*rand
end

% Vehicle trajectory information
vx = 25;        % longitudinal speed (m/s)
p0 = [0,2];     % initial position of vehicle (m,m)
h0 = 0;      % initial heading of vehicle (rad)
t_h = 14.7;         % time horizon to check (s)
% Create a vector of the trajectory information
x0 = [p0'; h0; vx];

% Vehicle dimensional information
vehicle.dr = 1.9;       % CG-front bumper distance (m)
vehicle.df = 3.0;        % CG-rear bumper distance (m)
vehicle.w = 1.8;        % vehicle width (m)


%%

for j = 1:size(barrellocations,1)
    angles(j,1) = 2*pi*rand;
end

% converting points to barrels for channeling device
barrelsArray2 = fcn_Patch_fillBarrelAtPoint(barrellocations,angles);

% detect collision and show visualization
[collFlags,collTime,collLoc,clearance,bodyCollLoc] = fcn_Patch_checkStraightCollisions(x0,vehicle,barrelsArray2,t_h,fig_num2);

