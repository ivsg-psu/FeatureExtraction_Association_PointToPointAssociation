function [collFlag,time,angle,location,clearance,bodyLoc] = fcn_Patch_checkCollisions(x0,vehicle,patchArray)
% fcn_Patch_checkCollisions
% Evaluates a circular vehicle trajectory against a series of patches to
% determine whether there will be collisions between the vehicle and the
% outline of any of the patch objects.
%
% ASSUMPTIONS:
%       1) The vehicle moves at constant speed along a circular trajectory.
%       2) The vehicle is represented by a bounding rectangle.
%       4) A positive trajectory radius indicates CCW travel. Negative
%       indicates CW.
%
% FORMAT:
%
%       [time,angle,location,clearance,bodyLoc] = fcn_Patch_checkCollisions(x0,vehicle,patchArray)
%
% INPUTS:
%
%      x0: a 6 x 1 vector containing the starting (x,y) coordinates of the
%           vehicle, the initial heading, the body slip angle, the
%           longitudinal vehicle speed, and the signed trajectory radius in
%           (m,m), radians, radians, m/s, and m.
%      vehicle: a structure containing the vehicle properties, which must
%           include fields df, dr, w for the vehicle CG-front bumper
%           distance, CG-rear bumper distance, and body width, in meters,
%           respectively.
%      patchArray: a structure array defining the objects with which the
%           vehicle could potentially collide
%
% OUTPUTS:
%
%      collFlag: an N x 1 vector of flags denoting whether there is a
%           collision with each of the objects
%      time: an N x 1 vector of collision times, where N is the number of
%           patch objects. Elements of the time vector will be set to Inf
%           if there is no overlap with the vehicle path.
%      angle: a N x 1
%      location: a N x 2 vector of collision locations, where N is the
%           number of patch objects and the columns are the x and y
%           coordinates of the collision. Elements of the location matrix
%           will be set to NaN if there is no overlap with the vehicle
%           path.
%       clearance: an N x 1 vector of minimum clearance distances between
%           the vehicle and the patch objects, where N is the number of
%           patch objects. Elements of the clearance vector will be set to
%           NaN if there is a collision with the object.
%       bodyLoc: an N x 2 vector of collision locations in vehicle body
%           fixed coordinates, where N is the number of patch objects and
%           the columns are the x and y coordinates of the collision.
%           Elements of the location matrix will be set to NaN if there is
%           no overlap with the vehicle path
%
% DEPENDENCIES:
%
%       fcn_geometry_findIntersectionLineSegmentWithCircle.m from the PSU
%       geometry library
%
% EXAMPLES:
%
%       See the script: script_test_fcn_Patch_checkCollisions.m for a full test
%       suite.
%
% This function was written by C. Beal
% Questions or comments? cbeal@bucknell.edu

% Revision history:
%     2022_02_17
%     -- wrote the code
%     2022_03_02
%     -- substantial debugging of odd cases

% TO DO
%   1) Does not work properly for negative radii
%   2) Search for other break cases
%       - collisions right after start points, hampered by determining
%       which of the double intersections to treat

flag_do_debug = 1; % Flag to plot the results for debugging
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
    if nargin ~= 3
        error('Incorrect number of input arguments')
    end
    
    % Check to see that there were any points provided
    if isempty(patchArray)
        warning('Empty array of objects, nothing to do.');
        return
    end
    
    % Check the trajectory input
    if size(x0,1) ~= 6
        error('Vehicle trajectory missing elements, cannot calculate collisions.')
    end
    
    % Check the vehicle structure input to make sure that the dimensions a,
    % b, and d are all supplied
    if ~all(isfield(vehicle,{'df','dr','w'}))
        error('One or more necessary vehicle dimensions missing. Check inputs.')
    end
    if isfield(vehicle,'df') && isempty(vehicle.df)
        error('CG-front axle distance empty.')
    end
    if isfield(vehicle,'dr') && isempty(vehicle.df)
        error('CG-rear axle distance empty.')
    end
    if isfield(vehicle,'w') && isempty(vehicle.df)
        error('Vehicle width empty.')
    end
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

% Break out some variables for easier referencing
p0 = x0(1:2);   % Initial location of the vehicle
h0 = x0(3);     % Initial heading of the vehicle
a0 = x0(4);     % Vehicle body slip angle
v0 = x0(5);     % Initial (and constant) vehicle speed
R = x0(6);      % Radius of the circular trajectory

% Enumerate some variables for easier code reading
LF = 1;
RF = 2;
RR = 3;
LR = 4;
IT = 5; % Inside tangent edge
X = 1;
Y = 2;

% Determine the number of patches to check
Npatches = length(patchArray);

% Calculate center point of vehicle trajectory circle
pc(X) = p0(X) + R*cos(h0+pi/2);
pc(Y) = p0(Y) + R*sin(h0+pi/2);

% Determine the unsigned bounding radii for all portions of
% the vehicle as well as the radii for the front left and front right corners
% Calculate the various pertinent radii and corner points of the vehicle
[radii,vehicleBB,radiiFlags] = fcn_Patch_CalcCircularTrajectoryGeometry(x0,vehicle);

% Parse out which radii are which
R_LF = radii(LF);
R_RF = radii(RF);
R_RR = radii(RR);
R_LR = radii(LR);
R_IT = radii(IT);
Rmin = radii(6);
Rmax = radii(7);

% Determine which edges of the vehicle are relevant to check. A side needs
% to be checked if the rear point is on the min or max radius
flag_check_left_side = 1; flag_check_left_tangent = 0;
flag_check_right_side = 1; flag_check_right_tangent = 0;
% if LF == radiiFlags(1) || LF == radiiFlags(2)
%     flag_check_left_side = 0;
% end
if (R > 0 && IT == radiiFlags(1))
    flag_check_left_tangent = 1;
    flag_check_left_side = 0;
end
% if RF == radiiFlags(1) || RF == radiiFlags(2)
%     flag_check_right_side = 0;
% end
if (R < 0 && IT == radiiFlags(1))
    flag_check_right_tangent = 1;
    flag_check_right_side = 0;
end
if flag_do_debug
    if flag_check_left_side && flag_check_left_tangent
        error('Bad combination of vehicle side checks (left side)')
    end
    if flag_check_right_side && flag_check_right_tangent
        error('Bad combination of vehicle side checks (right side)')
    end
end

% Pre-allocate the results with negative numbers (not viable results) for
% debugging purposes
time = -ones(Npatches,1);
angle = -ones(Npatches,1);
location = -ones(Npatches,2);
clearance = -ones(Npatches,1);
collFlag = -ones(Npatches,1);
bodyLoc = -ones(Npatches,2);

% Iterate over all of the patches in the patchArray input
for patchInd = 1:Npatches
    % Pull out the points into temporary vectors for code brevity
    xobst = patchArray(patchInd).pointsX;
    yobst = patchArray(patchInd).pointsY;
    
    % Determine the distance from the center of the vehicle trajectory for
    % all vertices on the object
    vertexRadii = sqrt((xobst - pc(1)).^2 + (yobst - pc(2)).^2);
    vertexAngles = atan2(yobst - pc(2), xobst - pc(1));
    
    % Pre-allocate matrices to store collision locations
    Nvertices = size(xobst,1);
    thetaValues = nan(Nvertices,7);
    
    %thetaLeftFrontCorner = nan(Nvertices,1);
    xyLeftFrontCorner = nan(Nvertices,2);
    bodyXYLeftFrontCorner = nan(Nvertices,2);
    %thetaRightFrontCorner = nan(Nvertices,1);
    xyRightFrontCorner = nan(Nvertices,2);
    bodyXYRightFrontCorner = nan(Nvertices,2);
    %thetaLeftRearCorner = nan(Nvertices,1);
    xyLeftRearCorner = nan(Nvertices,2);
    bodyXYLeftRearCorner = nan(Nvertices,2);
    %thetaRightFrontCorner = nan(Nvertices,1);
    xyRightRearCorner = nan(Nvertices,2);
    bodyXYRightRearCorner = nan(Nvertices,2);
    %thetaLeftSide = nan(Nvertices,1);
    xyLeftSide = nan(Nvertices,2);
    bodyXYLeftSide = nan(Nvertices,2);
    %thetaRightSide = nan(Nvertices,1);
    xyRightSide = nan(Nvertices,2);
    bodyXYRightSide = nan(Nvertices,2);
    %thetaFront = nan(Nvertices,1);
    xyFront = nan(Nvertices,2);
    bodyXYFront = nan(Nvertices,2);
    
    % Loop through each of the object vertices and associated edges,
    % checking for intersections with the vehicle vertices and edges
    for vertexInd = 1:Nvertices
        % Determine the index of the next vertex to define edges
        nextVertex = mod(vertexInd,Nvertices)+1;
        % Determine the endpoints of the object edge to be checked
        pa = [xobst(vertexInd) yobst(vertexInd)];
        pb = [xobst(nextVertex) yobst(nextVertex)];
        
        %%%%% Section 1: vehicle vertex to object edge checks %%%%%
        % If the left front corner of the vehicle is exposed
        if flag_check_left_side || flag_check_left_tangent
            [intersectionAngles,intersectionPoint] = ...
                fcn_geometry_findIntersectionLineSegmentWithCircle(pa,pb,pc,R_LF);
            if 0 < size(intersectionAngles,1)
                bodyXYLeftFrontCorner(vertexInd,:) = [vehicle.df vehicle.w/2];
                thetaLeftFrontStart = atan2(vehicleBB(LF,Y)-pc(Y),vehicleBB(LF,X)-pc(X));
                [~,minInd] = min(intersectionAngles);
                thetaValues(vertexInd,LF) = intersectionAngles(minInd) - thetaLeftFrontStart;
                xyLeftFrontCorner(vertexInd,:) = intersectionPoint(minInd,:);
            end
        end
        % If the right front corner of the vehicle is exposed
        if flag_check_right_side || flag_check_right_tangent
            [intersectionAngles,intersectionPoint] = ...
                fcn_geometry_findIntersectionLineSegmentWithCircle(pa,pb,pc,R_RF);
            if 0 < size(intersectionAngles,1)
                bodyXYRightFrontCorner(vertexInd,:) = [vehicle.df -vehicle.w/2];
                thetaRightRearStart = atan2(vehicleBB(RF,Y)-pc(Y),vehicleBB(RF,X)-pc(X));
                [~,minInd] = min(intersectionAngles);
                thetaValues(vertexInd,RF) = intersectionAngles(minInd) - thetaRightRearStart;
                xyRightFrontCorner(vertexInd,:) = intersectionPoint(minInd,:);
            end
        end
        % If the left rear corner of the vehicle is exposed
        if flag_check_left_side
            [intersectionAngles,intersectionPoint] = ...
                fcn_geometry_findIntersectionLineSegmentWithCircle(pa,pb,pc,R_LR);
            if 0 < size(intersectionAngles,1)
                bodyXYLeftRearCorner(vertexInd,:) = [-vehicle.dr vehicle.w/2];
                thetaLeftRearStart = atan2(vehicleBB(LR,Y)-pc(Y),vehicleBB(LR,X)-pc(X));
                [~,minInd] = min(intersectionAngles);
                thetaValues(vertexInd,LR) = intersectionAngles(minInd) - thetaLeftRearStart;
                xyLeftRearCorner(vertexInd,:) = intersectionPoint(minInd,:);
            end
        end
        % If the right rear corner of the vehicle is exposed
        if flag_check_right_side
            [intersectionAngles,intersectionPoint] = ...
                fcn_geometry_findIntersectionLineSegmentWithCircle(pa,pb,pc,R_RR);
            if 0 < size(intersectionAngles,1)
                bodyXYRightRearCorner(vertexInd,:) = [-vehicle.dr -vehicle.w/2];
                thetaRightRearStart = atan2(vehicleBB(RR,Y)-pc(Y),vehicleBB(RR,X)-pc(X));
                [~,minInd] = min(intersectionAngles);
                thetaValues(vertexInd,RR) = intersectionAngles(minInd) - thetaRightRearStart;
                xyRightRearCorner(vertexInd,:) = intersectionPoint(minInd,:);
            end
        end
        
        %%%%% Section 2: object vertex to vehicle edge checks %%%%%
        
        % First check to see if the vertex falls into the vehicle
        % trajectory at all
        if vertexRadii(vertexInd) > Rmin && vertexRadii(vertexInd) < Rmax
            % Next, see if the left side of the vehicle is exposed
            if flag_check_left_side || flag_check_left_tangent
                if vertexRadii(vertexInd) < R_LF
                    % vertex has an intersection with the left side of the
                    % vehicle
                    if flag_check_left_tangent
                        % vertex has an intersection ahead of the tangent
                        % point
                        [thetaOffset,~] = ...
                            fcn_geometry_findIntersectionLineSegmentWithCircle(vehicleBB(IT,:),vehicleBB(LF,:),pc,vertexRadii(vertexInd));
                    else
                        % vertex has an intersection ahead of the LR vertex
                        [thetaOffset,~] = ...
                            fcn_geometry_findIntersectionLineSegmentWithCircle(vehicleBB(LR,:),vehicleBB(LF,:),pc,vertexRadii(vertexInd));
                    end
                    if ~isempty(thetaOffset)
                        thetaValues(vertexInd,5) = vertexAngles(vertexInd)-thetaOffset;
                        xyLeftSide(vertexInd,:) = [xobst(vertexInd),yobst(vertexInd)];
                    end
                    % CEB: Need to determine the body xy coordinates for
                    % these cases
                end
            end
            % Next, see if the left side of the vehicle is exposed
            if flag_check_right_side || flag_check_right_tangent
                if vertexRadii(vertexInd) > R_RF
                    % vertex has an intersection with the right side of the
                    % vehicle
                    if flag_check_right_tangent
                        % vertex has an intersection ahead of the tangent
                        % point
                        [thetaOffset,~] = ...
                            fcn_geometry_findIntersectionLineSegmentWithCircle(vehicleBB(IT,:),vehicleBB(RF,:),pc,vertexRadii(vertexInd));
                    else
                        % vertex has an intersection ahead of the RR vertex
                        [thetaOffset,~] = ...
                            fcn_geometry_findIntersectionLineSegmentWithCircle(vehicleBB(RR,:),vehicleBB(RF,:),pc,vertexRadii(vertexInd));
                    end
                    if ~isempty(thetaOffset)
                        thetaValues(vertexInd,6) = vertexAngles(vertexInd)-thetaOffset;
                        xyRightSide(vertexInd,:) = [xobst(vertexInd),yobst(vertexInd)];
                    end
                    % CEB: Need to determine the body xy coordinates for
                    % these cases
                end
            end
            % Finally, check the front of the vehicle
            if vertexRadii(vertexInd) >= R_LF && vertexRadii(vertexInd) <= R_RF
                % vertex has an intersection with the front of the vehicle
                [thetaOffset,~] = fcn_geometry_findIntersectionLineSegmentWithCircle(vehicleBB(LF,:),vehicleBB(RF,:),pc,vertexRadii(vertexInd));
                thetaValues(vertexInd,7) = vertexAngles(vertexInd)-thetaOffset;
                xyFront(vertexInd,:) = [xobst(vertexInd),yobst(vertexInd)];
                % CEB: Need to determine the body xy coordinates for this
                % case
            end
        end
    end
    
    %%%%% Section 3: first collision selection, or clearance calc %%%%%
    
    % Check for a case where no intersections were calculated
    if all(isnan(thetaValues),'all')
        collFlag(patchInd) = 0;
        
        % Determine the nearest two vertices
        [nearestOuterClearance,nearestOuters] = mink(min(inf(Nvertices,1),vertexRadii - abs(Rmax)),2);
        [nearestInnerClearance,nearestInners] = mink(min(inf(Nvertices,1),abs(Rmin) - vertexRadii),2);
        if max(0,nearestOuterClearance) > max(nearestInnerClearance)
            pa = [xobst(nearestOuters(1)) yobst(nearestOuters(1))];
            pb = [xobst(nearestOuters(2)) yobst(nearestOuters(2))];
            theta_offset = sign(R)*asin(vehicle.dr/abs(R));
            % Compute the alpha associated with the min radius
            alpha = -((pa(1)-pc(1))*(pb(1)-pa(1)) + (pa(2)-pc(2))*(pb(2)-pa(2)))/((pb(1)-pa(1))^2 + (pb(2)-pa(2))^2);
            % Check to see if the minimum clearance is along the edge between
            % the nearest two vertices
            if alpha > 0 && alpha < 1
                % If so, calculate the clearance distance from the computed
                % point to the
                clearance(patchInd) = norm(pa + alpha*(pb-pa) - pc) - R;
                location(patchInd,:) = pa + alpha*(pb-pa);
            else
                clearance(patchInd) = norm(pa - pc) - R;
                location(patchInd,:) = pa;
            end
        else
            pa = [xobst(nearestInners(1)) yobst(nearestInners(1))];
            theta_offset = 0;
            clearance(patchInd) = abs(norm(pa - pc) - R);
            location(patchInd,:) = pa;
        end
        angle(patchInd) = atan2(location(patchInd,2)-pc(2),location(patchInd,1)-pc(1)) + theta_offset;
        if R >= 0
            time(patchInd) = rerangeAngles(angle(patchInd) - h0 + pi/2)*abs(R)/v0;
        else
            time(patchInd) = rerangeAngles(h0 + pi/2 - angle(patchInd))*abs(R)/v0;
        end
        % No collision, so no body collision location
        bodyLoc(patchInd,:) = [NaN NaN];
        
        % Now, find the minimum angular location for the patch object in the
        % vehicle travel direction (indicated by the sign of R)
    else
        if R >= 0
            collFlag(patchInd) = 1;
            % Find the index with the minimum angle in each dimension. Then,
            % that identifies the single case which yields the nearest
            % collision. Pull the data appropriately.
            thetaValues = rerangeAngles(thetaValues);
            [minVals,idx] = nanmin(thetaValues);
            [~,idy] = nanmin(minVals);
            idx = idx(idy);
        else
            collFlag(patchInd) = 1;
            % Find the index with the minimum angle in each dimension. Then,
            % that identifies the single case which yields the nearest
            % collision. Pull the data appropriately.
            thetaValues = rerangeAngles(2*pi-thetaValues);
            [minVals,idx] = nanmin(thetaValues);
            [~,idy] = nanmin(minVals);
            idx = idx(idy);
        end
        
        % Use the indices of the smallest angle value to pull the correct
        % collision location data for output
        angle(patchInd) = thetaValues(idx,idy);
        
        switch(idy)
            case LF
                location(patchInd,:) = xyLeftFrontCorner(idx,:);
                bodyLoc(patchInd,:) = bodyXYLeftFrontCorner(idx,:);
            case RF
                location(patchInd,:) = xyRightFrontCorner(idx,:);
                bodyLoc(patchInd,:) = bodyXYRightFrontCorner(idx,:);
            case RR
                location(patchInd,:) = xyRightRearCorner(idx,:);
                bodyLoc(patchInd,:) = bodyXYRightRearCorner(idx,:);
            case LR
                location(patchInd,:) = xyLeftRearCorner(idx,:);
                bodyLoc(patchInd,:) = bodyXYLeftRearCorner(idx,:);
            case 5
                location(patchInd,:) = xyLeftSide(idx,:);
                bodyLoc(patchInd,:) = bodyXYLeftSide(idx,:);
            case 6
                location(patchInd,:) = xyRightSide(idx,1);
                bodyLoc(patchInd,:) = bodyXYRightSide(idx,:);
            case 7
                location(patchInd,:) = xyFront(idx,:);
                bodyLoc(patchInd,:) = bodyXYFront(idx,:);
        end
        
        % With the location set, determine the time required to reach the
        % location
        time(patchInd) = (angle(patchInd) - h0 + pi/2)*abs(R)/v0;
    end
end
end

% Function to re-range angles between 0 and 2*pi
function outAngles = rerangeAngles(inAngles)
outAngles = inAngles;
while any(outAngles > 2*pi,'all')
    outAngles(outAngles > 2*pi) = outAngles(outAngles > 2*pi) - 2*pi;
end
while any(outAngles < 0,'all')
    outAngles(outAngles < 0) = outAngles(outAngles < 0) + 2*pi;
end
end