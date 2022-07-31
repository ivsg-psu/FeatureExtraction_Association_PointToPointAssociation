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

bottomLaneMarker = [0 0.1; 367 0.1];         % Bottom lane marker
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
% 
% % Second points for  row of channeling device
% barrels(19,1) = 0; barrels(19,2) = 11.1;
% for i = 19:35
%     barrels(i+1,1) = barrels(i) + 18.5;
%     barrels(i+1,2) = 11.1;
% end
% 
% % left diagonal points of channeling device
% for i = 36:40
%     barrels(i+1,1) = barrels(19,1) + (i-35)*8.3;
%     barrels(i+1,2) = barrels(19,2) - (i-35)*0.55;
% end
% 
% % right diagonal points of channeling device
% for i = 41:45
%     barrels(i+1,1) = barrels(36,1) + (i-40)*8.3;
%     barrels(i+1,2) = barrels(36,2) - (i-40)*0.55;
% end

angles = zeros([46 2]);
% randomized orientation of the barrels
for j = 1:46
    angles(j,1) = 2*pi*rand;
end

% converting points to barrels for channeling device
barrelsArray = fcn_Patch_fillBarrelAtPoint(barrels,angles);

% detect collision and show visualization
[collFlags,collTime,collLoc,clearance,bodyCollLoc] = fcn_Patch_checkStraightCollisions(x0,vehicle,barrelsArray,t_h,fig_num);

%% Example 2
%
clear, close

fig_num2 = 2;
figure(fig_num2)
hold on

w = 3.3528;  % Lane width in meters
l = 353.368; % Total lane length in meters

blockedLaneMarker = [0 11.1; 353.568 11.1];       % Lane line to represent blocked lane
plot(blockedLaneMarker(:,1),blockedLaneMarker(:,2),'k')
axis([0 367 -150 150])

t_w = (0.1+2*w); % top lane y position
topLaneMarker = [0 t_w; l t_w];           % Top lane marker
plot(topLaneMarker(:,1),topLaneMarker(:,2),'k--')

m_w = (0.1+w); % middle lane y position
middleLaneMarker = [0 (0.1+m_w); l m_w];        % Middle lane marker
plot(middleLaneMarker(:,1),middleLaneMarker(:,2),'k--')

bottomLaneMarker = [0 0.1; l 0.1];         % Bottom lane marker
plot(bottomLaneMarker(:,1),bottomLaneMarker(:,2),'k')

bc_w = (0.1-w); % bottom closed lane y position
bottomClosedLaneMarker = [0 bc_w; l bc_w];         % Bottom closed lane marker
plot(bottomClosedLaneMarker(:,1),bottomClosedLaneMarker(:,2),'k')

% points for diagonal of channeling device
barrels2(1,1) = 15.24; barrels2(1,2) = m_w;
for i = 1:5  
    barrels2(i+1,1) = barrels2(i,1) + 3.048;
    barrels2(i+1,2) = barrels2(i,2) + w/5;
end

% points for row of channeling device
for i = 6:14
    barrels2(i+1,1) = barrels2(i,1) + 18.288;
    barrels2(i+1,2) = t_w;
end


% points for diagonal of channeling device
barrels2(16,1) = 210.312; barrels2(16,2) = t_w;
for i = 16:20
    barrels2(i+1,1) = barrels2(i,1) + 4.572;
    barrels2(i+1,2) = barrels2(i,2) - w/5;
end

% points for row of channeling device
barrels2(22,1) = barrels2(15,1) + 36.576; barrels2(22,2) = t_w;
for i = 22:24
    barrels2(i+1,1) = barrels2(i,1) + 18.288;
    barrels2(i+1,2) = t_w;
end

% points for diagonal of channeling device
barrels2(26,1) = 298.704; barrels2(26,2) = t_w;
for i = 26: 31
    barrels2(i+1,1) = barrels2(i,1) + 9.144;
    barrels2(i+1,2) = barrels2(i,2) - w/6;
end

barrels2(33,1) = 118.872; barrels2(33,2) = m_w;
for i = 33: 37
    barrels2(i+1,1) = barrels2(i,1) + 18.288;
    barrels2(i+1,2) = m_w;
end

for i = 38: 42
    barrels2(i+1,1) = barrels2(i,1) + 15.24/5;
    barrels2(i+1,2) = barrels2(i,2) - w/5;
end

for i = 43:52
    barrels2(i+1,1) = barrels2(33,1) - 45.72*(i-42)/11;
    barrels2(i+1,2) = barrels2(33,2) - 2*w*(i-42)/11;
end

for i = 53:58
    barrels2(i+1,1) = barrels2(21,1) + 115.824*(i-52)/6;
    barrels2(i+1,2) = m_w;
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

angles2 = zeros([size(barrels2,1) 2]);
for j = 1:size(barrels2,1)
    angles2(j,1) = 2*pi*rand;
end

% converting points to barrels for channeling device
barrelsArray2 = fcn_Patch_fillBarrelAtPoint(barrels2,angles2);

fcn_Patch_plotPatch(barrelsArray2,2)
axis square

% % detect collision and show visualization
% [collFlags,collTime,collLoc,clearance,bodyCollLoc] = fcn_Patch_checkStraightCollisions(x0,vehicle,barrelsArray2,t_h,fig_num2);

%% Example 3
% Same scenario as above with radial noise

fig_num3 = 3;
figure(fig_num3)
hold on

blockedLaneMarker = [0 11.1; 353.568 11.1];       % Lane line to represent blocked lane
plot(blockedLaneMarker(:,1),blockedLaneMarker(:,2),'k')
axis([0 367 -150 150])

t_w = (0.1+2*w); % top lane y position
topLaneMarker = [0 t_w; l t_w];           % Top lane marker
plot(topLaneMarker(:,1),topLaneMarker(:,2),'k--')

m_w = (0.1+w); % middle lane y position
middleLaneMarker = [0 (0.1+m_w); l m_w];        % Middle lane marker
plot(middleLaneMarker(:,1),middleLaneMarker(:,2),'k--')

bottomLaneMarker = [0 0.1; l 0.1];         % Bottom lane marker
plot(bottomLaneMarker(:,1),bottomLaneMarker(:,2),'k')

bc_w = (0.1-w); % bottom closed lane y position
bottomClosedLaneMarker = [0 bc_w; l bc_w];         % Bottom closed lane marker
plot(bottomClosedLaneMarker(:,1),bottomClosedLaneMarker(:,2),'k')

dataset = {barrels2};
newbarrels2 = fcn_Points_addRadialNoise(dataset,0.1);

new_angles2 = zeros([size(barrels2,1) 2]);
for j = 1:size(barrels2,1)
    new_angles2(j,1) = 2*pi*rand;
end

barrelsArray2 = fcn_Patch_fillBarrelAtPoint(newbarrels2{1},new_angles2);

fcn_Patch_plotPatch(barrelsArray2,3)
axis square

%% Example 3 continued
% performing point association
[pairedXYdata, numMatches, nonMatchesA, nonMatchesB]  = fcn_Points_pairXYdata(newbarrels2{1},barrels2,1);

figure(4);
clf
hold on
grid on

for i = 1:numMatches
    plot([pairedXYdata(i,1) pairedXYdata(i,3)],[pairedXYdata(i,2) pairedXYdata(i,4)],'o','MarkerSize',10,'MarkerFaceColor',[0 1 0])
end
quiver(pairedXYdata(:,1),pairedXYdata(:,2),pairedXYdata(:,3)-pairedXYdata(:,1),pairedXYdata(:,4)-pairedXYdata(:,2),0,'k','linewidth',1)
% Also plot the non-matches in data set A with circles in various colors
for i = 1+numMatches:nonMatchesA+numMatches
    plot(pairedXYdata(i,1),pairedXYdata(i,2),'o','MarkerSize',10,'MarkerFaceColor',[1 0 0])
end
% Lastly, plot the non-matches in data set B with squares in various colors
for i = 1+numMatches+nonMatchesA:nonMatchesA+numMatches+nonMatchesA
    plot(pairedXYdata(i,3),pairedXYdata(i,4),'o','MarkerSize',10,'MarkerFaceColor',[0 0 1])
end