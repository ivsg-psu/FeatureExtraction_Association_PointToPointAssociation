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
%       [collFlag,time,angle,location,clearance,bodyLoc] = fcn_Patch_checkCollisions(x0,vehicle,patchArray)
%
% INPUTS:
%
%      x0: a 6 x 1 vector containing the starting (x,y) coordinates of the
%           vehicle, the initial course, the body slip angle, the
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
%       (https://github.com/ivsg-psu/PathPlanning_GeomTools_GeomClassLibrary.git)
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
%     2022_03_15
%     -- fixed intersection detection and near-miss position calculation,
%     added additional comments

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
    if nargin ~= 3
        error('Incorrect number of input arguments')
    end
    
    % Check to see that there were any obstacle patch structures
    if isempty(patchArray)
        warning('Empty array of objects, nothing to do.');
        return
    end
    
    % Check the trajectory input
    if size(x0,1) ~= 6
        error('Vehicle trajectory missing elements, cannot calculate collisions.')
    end
    
    % Check the vehicle structure input to make sure that the dimensions w,
    % df, and dr are all supplied
    if ~all(isfield(vehicle,{'df','dr','w'}))
        error('One or more necessary vehicle dimensions missing. Check inputs.')
    end
    if isfield(vehicle,'df') && isempty(vehicle.df)
        error('CG-front bumper distance empty.')
    end
    if isfield(vehicle,'dr') && isempty(vehicle.df)
        error('CG-rear bumper distance empty.')
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
CoG0 = x0(3);     % Initial heading of the vehicle
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
pc(X) = p0(X) + R*cos(CoG0+pi/2);
pc(Y) = p0(Y) + R*sin(CoG0+pi/2);

% Determine the unsigned bounding radii for all portions of
% the vehicle as well as the radii for the front left and front right corners
% Calculate the various pertinent radii and corner points of the vehicle
[radii,vehicleBB,radiiFlags] = fcn_Patch_CalcCircularTrajectoryGeometry(x0,vehicle);

% Parse out which radii are which from the trajectory calculation
R_LF = radii(LF);
R_RF = radii(RF);
R_RR = radii(RR);
R_LR = radii(LR);
R_IT = radii(IT);
Rmin = radii(6);
Rmax = radii(7);

% Determine whether the vehicle has a left or right tangent to the minimum
% radius. Set the flags to indicate that there is not one (minimum is on a
% vehicle corner) initially.
flag_check_left_tangent = 0;
flag_check_right_tangent = 0;

% Check and see if the minimum radius is tangent to left side, based on the
% sign of the radius and the flag from the circular geometry calculations
if (R > 0 && IT == radiiFlags(1))
    flag_check_left_tangent = 1;
end
% Check and see if the minimum radius is tangent to right side, based on
% the sign of the radius and the flag from the circular geometry calculations
if (R < 0 && IT == radiiFlags(1))
    flag_check_right_tangent = 1;
end

% Pre-allocate the results with Nans (not viable results) for
% debugging purposes
time = nan(Npatches,1);
angle = nan(Npatches,1);
location = nan(Npatches,2);
clearance = nan(Npatches,1);
collFlag = nan(Npatches,1);
bodyLoc = nan(Npatches,2);

thetaLeftFrontStart = atan2(vehicleBB(LF,Y)-pc(Y),vehicleBB(LF,X)-pc(X));
thetaRightFrontStart = atan2(vehicleBB(RF,Y)-pc(Y),vehicleBB(RF,X)-pc(X));
thetaLeftRearStart = atan2(vehicleBB(LR,Y)-pc(Y),vehicleBB(LR,X)-pc(X));
thetaRightRearStart = atan2(vehicleBB(RR,Y)-pc(Y),vehicleBB(RR,X)-pc(X));
thetaITStart = atan2(vehicleBB(5,Y)-pc(Y),vehicleBB(5,X)-pc(X));
        
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
    % There are 14 columns to store two potential intersections for each
    % possible check: each of the four vehicle vertices to an object edge
    % and an object vertex to each of the potentially exposed vehicle sides
    % (left, right, front). This could be expanded to N vertices on the
    % vehicle as well, for more complex vehicle shapes.
    thetaValues = nan(Nvertices,14);
    
    % Allocate space to store the intersection locations on the vehicle
    % body and in inertial space. The corners of the vehicle are known, so
    % set the values of those right away. Note that this changes the code
    % slightly when the results are processed for output later.
    bodyXYLeftFrontCorner = [vehicle.df vehicle.w/2];
    bodyXYRightFrontCorner = [vehicle.df -vehicle.w/2];
    bodyXYRightRearCorner = [-vehicle.dr -vehicle.w/2];
    bodyXYLeftRearCorner = [-vehicle.dr vehicle.w/2];
    % Allocate the rest with NaN values that can be ignored if not filled
    xyLeftFrontCorner = nan(Nvertices,4);
    xyRightFrontCorner = nan(Nvertices,4);
    xyRightRearCorner = nan(Nvertices,4);
    xyLeftRearCorner = nan(Nvertices,4);
    xyLeftSide = nan(Nvertices,4);
    bodyXYLeftSide = nan(Nvertices,4);
    xyRightSide = nan(Nvertices,4);
    bodyXYRightSide = nan(Nvertices,4);
    xyFront = nan(Nvertices,4);
    bodyXYFront = nan(Nvertices,4);
    
    % Loop through each of the object vertices and associated edges,
    % checking for intersections with the vehicle vertices and edges
    for vertexInd = 1:Nvertices
        % Determine the index of the next vertex to define edges
        nextVertex = mod(vertexInd,Nvertices)+1;
        % Determine the endpoints of the object edge to be checked
        pa = [xobst(vertexInd) yobst(vertexInd)];
        pb = [xobst(nextVertex) yobst(nextVertex)];
        
        %%%%% Section 1: vehicle vertex to object edge checks %%%%%
        
        % Check the left front corner of the vehicle for intersections with
        % the currently selected obstacle edge (between pa and pb)
        [intersectionAngles,intersectionPoint] = ...
            fcn_geometry_findIntersectionLineSegmentWithCircle(pa,pb,pc,R_LF);
        % Iterate through any intersections found
        for intIndex = 1:size(intersectionAngles,2)
            % Compute the angle of the intersection and subtract the angle
            % of the LF point at the start of the vehicle trajectory to get
            % the traveled angular distance to the collision point
            thetaValues(vertexInd,2*(LF-1)+intIndex) = intersectionAngles(intIndex) - thetaLeftFrontStart;
            % Store the calculated intersection point
            xyLeftFrontCorner(vertexInd,(2*intIndex-1):(2*intIndex)) = intersectionPoint(intIndex,:);
            % Plot, if flagged
            if flag_do_debug
                plot(xyLeftFrontCorner(vertexInd,2*intIndex-1),xyLeftFrontCorner(vertexInd,2*intIndex),'rd')
            end
        end
        % Check the right front corner of the vehicle for intersections with
        % the currently selected obstacle edge (between pa and pb)
        [intersectionAngles,intersectionPoint] = ...
            fcn_geometry_findIntersectionLineSegmentWithCircle(pa,pb,pc,R_RF);
        % Iterate through any intersections found
        for intIndex = 1:size(intersectionAngles,2)
            % Compute the angle of the intersection and subtract the angle
            % of the RF point at the start of the vehicle trajectory to get
            % the traveled angular distance to the collision point
            thetaValues(vertexInd,2*(RF-1)+intIndex) = intersectionAngles(intIndex) - thetaRightFrontStart;
            % Store the calculated intersection point
            xyRightFrontCorner(vertexInd,(2*intIndex-1):(2*intIndex)) = intersectionPoint(intIndex,:);
            % Plot, if flagged
            if flag_do_debug
                plot(xyRightFrontCorner(vertexInd,2*intIndex-1),xyRightFrontCorner(vertexInd,2*intIndex),'bd')
            end
        end
        % If the left side is tangent to the minimum radius, check the
        % tangent point for intersections with the currently selected
        % obstacle edge (between pa and pb)
        if flag_check_left_tangent
            [intersectionAngles,intersectionPoint] = ...
                fcn_geometry_findIntersectionLineSegmentWithCircle(pa,pb,pc,R_IT);
            theta_offset = thetaITStart;
        % If the left side is NOT tangent to the minimum radius, check the
        % LR point for intersections with the currently selected obstacle
        % edge (between pa and pb). (Intersections could happen behind the
        % front point in a "sideswipe" scenario.)
        else
            [intersectionAngles,intersectionPoint] = ...
                fcn_geometry_findIntersectionLineSegmentWithCircle(pa,pb,pc,R_LR);
            theta_offset = thetaLeftRearStart;
        end
        % Iterate through any intersections found
        for intIndex = 1:size(intersectionAngles,2)
            % Compute the angle of the intersection and subtract the angle
            % of the LR or tangent point at the start of the vehicle
            % trajectory to get the traveled angular distance to the
            % collision point
            thetaValues(vertexInd,2*(LR-1)+intIndex) = intersectionAngles(intIndex) - theta_offset;
            % Store the calculated intersection point
            xyLeftRearCorner(vertexInd,(2*intIndex-1):(2*intIndex)) = intersectionPoint(intIndex,:);
            % Plot, if flagged
            if flag_do_debug
                plot(xyLeftRearCorner(vertexInd,2*intIndex-1),xyLeftRearCorner(vertexInd,2*intIndex),'cd')
            end
        end
        % If the right side is tangent to the minimum radius, check the
        % tangent point for intersections with the currently selected
        % obstacle edge (between pa and pb)
        if flag_check_right_tangent
            [intersectionAngles,intersectionPoint] = ...
                fcn_geometry_findIntersectionLineSegmentWithCircle(pa,pb,pc,R_IT);
            theta_offset = thetaITStart;
        % If the right side is NOT tangent to the minimum radius, check the
        % RR point for intersections with the currently selected obstacle
        % edge (between pa and pb). (Intersections could happen behind the
        % front point in a "sideswipe" scenario.)
        else
            [intersectionAngles,intersectionPoint] = ...
                fcn_geometry_findIntersectionLineSegmentWithCircle(pa,pb,pc,R_RR);
            theta_offset = thetaRightRearStart;
        end
        % Iterate through any intersections found
        for intIndex = 1:size(intersectionAngles,2)
            % Compute the angle of the intersection and subtract the angle
            % of the RR or tangent point at the start of the vehicle
            % trajectory to get the traveled angular distance to the
            % collision point
            thetaValues(vertexInd,2*(RR-1)+intIndex) = intersectionAngles(intIndex) - theta_offset;
            % Store the calculated intersection point
            xyRightRearCorner(vertexInd,(2*intIndex-1):(2*intIndex)) = intersectionPoint(intIndex,:);
            % Plot, if flagged
            if flag_do_debug
                plot(xyRightRearCorner(vertexInd,2*intIndex-1),xyRightRearCorner(vertexInd,2*intIndex),'md')
            end
        end
        
        
        %%%%% Section 2: object vertex to vehicle edge checks %%%%%
        
        % First check to see if the vertex falls into the vehicle
        % trajectory at all
        if vertexRadii(vertexInd) > Rmin && vertexRadii(vertexInd) < Rmax
            % First, check the front of the vehicle, which will have a
            % collision if the obstacle vertex lies between the radii of
            % the LF and RF points. But need to check both cases since
            % the ordering of the radii magnitude will change depending on
            % the direction of travel.
            if (vertexRadii(vertexInd) >= R_LF && vertexRadii(vertexInd) <= R_RF) || (vertexRadii(vertexInd) <= R_LF && vertexRadii(vertexInd) >= R_RF)
                % If the vertex has an intersection with the front of the
                % vehicle, calculate the intersections
                [intersectionAngles,~] = ...
                    fcn_geometry_findIntersectionLineSegmentWithCircle(vehicleBB(LF,:),vehicleBB(RF,:),pc,vertexRadii(vertexInd));
                % Iterate through any intersections found
                for intIndex = 1:size(intersectionAngles,2)
                    % Calculate the angle of the vertex minus the angle on
                    % the at the starting position where the radius of the
                    % vertex hits the front of the vehicle. This yields the
                    % angle of the CG at the collision. Store any results
                    % in the 12th and 13th columns of the angles matrix.
                    thetaValues(vertexInd,12+intIndex) = vertexAngles(vertexInd)-intersectionAngles(intIndex);
                    % Store the calculated intersection point(s) in the 1st
                    % and 2nd or 3rd and 4th columns of the collision
                    % location matrix
                    xyFront(vertexInd,(2*intIndex-1):(2*intIndex)) = [xobst(vertexInd),yobst(vertexInd)];
                    % Plot, if flagged
                    if flag_do_debug
                        plot(xyFront(vertexInd,2*intIndex-1),xyFront(vertexInd,2*intIndex),'gs')
                    end
                    % CEB: Need to determine the body xy coordinates for this case
                end
            end
            % Next, calculate any itersections with the left side of the
            % vehicle. Start by determining how much of the left side of
            % the vehicle to check. For an inner side of the vehicle, only
            % the portion ahead of the tangent point needs to be checked,
            % if a tangent exists.
            if flag_check_left_tangent
                % A tangent condition exists, so calculate the intersection
                % angle relative to the starting point of the vehicle using
                % the portion of the left side ahead of the tangent point
                [intersectionAngles,~] = ...
                    fcn_geometry_findIntersectionLineSegmentWithCircle(vehicleBB(IT,:),vehicleBB(LF,:),pc,vertexRadii(vertexInd));
            else
                % No tangent condition exists, so calculate the
                % intersection angle relative to the starting point of the
                % vehicle using the entire left side of the vehicle
                [intersectionAngles,~] = ...
                    fcn_geometry_findIntersectionLineSegmentWithCircle(vehicleBB(LR,:),vehicleBB(LF,:),pc,vertexRadii(vertexInd));
            end
            % Iterate through any intersections found
            for intIndex = 1:size(intersectionAngles,2)
                % Calculate the angle of the vertex minus the angle on the
                % at the starting position where the radius of the vertex
                % hits the front of the vehicle. This yields the angle of
                % the CG at the collision. Store any results in the 9th
                % and 10th columns of the angles matrix.
                thetaValues(vertexInd,8+intIndex) = vertexAngles(vertexInd)-intersectionAngles(intIndex);
                % Store the calculated intersection point(s) in the 1st and
                % 2nd or 3rd and 4th columns of the collision location
                % matrix
                xyLeftSide(vertexInd,(2*intIndex-1):(2*intIndex)) = [xobst(vertexInd),yobst(vertexInd)];
                % Plot, if flagged
                if flag_do_debug
                    plot(xyLeftSide(vertexInd,2*intIndex-1),xyLeftSide(vertexInd,2*intIndex),'rs')
                end
            end
            % CEB: Need to determine the body xy coordinates for these cases
            
            % Next, calculate any itersections with the right side of the
            % vehicle. Start by determining how much of the right side of
            % the vehicle to check. For an inner side of the vehicle, only
            % the portion ahead of the tangent point needs to be checked,
            % if a tangent exists.
            if flag_check_right_tangent
                % A tangent condition exists, so calculate the intersection
                % angle relative to the starting point of the vehicle using
                % the portion of the right side ahead of the tangent point
                [intersectionAngles,~] = ...
                    fcn_geometry_findIntersectionLineSegmentWithCircle(vehicleBB(IT,:),vehicleBB(RF,:),pc,vertexRadii(vertexInd));
            else
                % No tangent condition exists, so calculate the
                % intersection angle relative to the starting point of the
                % vehicle using the entire right side of the vehicle
                [intersectionAngles,~] = ...
                    fcn_geometry_findIntersectionLineSegmentWithCircle(vehicleBB(RR,:),vehicleBB(RF,:),pc,vertexRadii(vertexInd));
            end
            % Iterate through any intersections found
            for intIndex = 1:size(intersectionAngles,2)
                % Calculate the angle of the vertex minus the angle on the
                % at the starting position where the radius of the vertex
                % hits the front of the vehicle. This yields the angle of
                % the CG at the collision. Store any results in the 11th
                % and 12th columns of the angles matrix.
                thetaValues(vertexInd,10+intIndex) = vertexAngles(vertexInd)-intersectionAngles(intIndex);
                % Store the calculated intersection point(s) in the 1st and
                % 2nd or 3rd and 4th columns of the collision location
                % matrix
                xyRightSide(vertexInd,(2*intIndex-1):(2*intIndex)) = [xobst(vertexInd),yobst(vertexInd)];
                % Plot, if flagged
                if flag_do_debug
                    plot(xyRightSide(vertexInd,2*intIndex-1),xyRightSide(vertexInd,2*intIndex),'bs')
                end
            end
            % CEB: Need to determine the body xy coordinates for these cases
        end
    end
    
    %%%%% Section 3: first collision selection, or clearance calc %%%%%
    
    % Check for a case where no intersections were calculated
    if all(isnan(thetaValues),'all')
        collFlag(patchInd) = 0;
        
        % Determine the nearest two vertices on both inside and outside of
        % the vehicle trajectory
        [~,nearestOuters] = mink(min(inf(Nvertices,1),vertexRadii - abs(Rmax)),2);
        [nearestInnerClearance,nearestInner] = mink(min(inf(Nvertices,1),abs(Rmin) - vertexRadii),1);
        
        % For the outer points, there could be a closer point along the
        % edge of the obstacle (only on the outside due to the convexity of
        % the outer trajectory). Determine the closest point on the line to
        % the center of the circle.
        
        % Grab the two closest outer points
        pa = [xobst(nearestOuters(1)) yobst(nearestOuters(1))];
        pb = [xobst(nearestOuters(2)) yobst(nearestOuters(2))];
        
        % Compute the alpha associated with the min radius
        alpha = -((pa(1)-pc(1))*(pb(1)-pa(1)) + (pa(2)-pc(2))*(pb(2)-pa(2)))/((pb(1)-pa(1))^2 + (pb(2)-pa(2))^2);
        % Check to see if the minimum clearance is along the edge between
        % the nearest two vertices
        if alpha > 0 && alpha < 1
            pmin = pa + alpha*(pb-pa);
            nearestOuterClearance = norm(pmin - pc) - R;
        else
            nearestOuterClearance = norm(pa - pc) - R;
        end
           
        % Check to see whether the closest point outside the trajectory is
        % closer than the closest point inside the trajectory
        if nearestOuterClearance < nearestInnerClearance
            
            % If the closest point is outside, the closest point on the
            % vehicle will be the corner on the max radius
            maxCorner = radiiFlags(2);
            theta_offset = atan2(vehicleBB(maxCorner,Y)-pc(Y),vehicleBB(maxCorner,X)-pc(X));
        
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
            minCorner = radiiFlags(1);
            theta_offset = atan2(vehicleBB(minCorner,Y)-pc(Y),vehicleBB(minCorner,X)-pc(X));
            clearance(patchInd) = abs(norm(pa - pc) - R);
            location(patchInd,:) = pa;
        end
        angle(patchInd) = atan2(location(patchInd,2)-pc(2),location(patchInd,1)-pc(1)) - theta_offset;
        if R >= 0
            time(patchInd) = rerangeAngles(angle(patchInd) - CoG0 + pi/2)*abs(R)/v0;
        else
            time(patchInd) = rerangeAngles(CoG0 + pi/2 - angle(patchInd))*abs(R)/v0;
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
            case 1
                location(patchInd,:) = xyLeftFrontCorner(idx,1:2);
                bodyLoc(patchInd,:) = bodyXYLeftFrontCorner(1,1:2);
            case 2
                location(patchInd,:) = xyLeftFrontCorner(idx,3:4);
                bodyLoc(patchInd,:) = bodyXYLeftFrontCorner(1,1:2);
            case 3
                location(patchInd,:) = xyRightFrontCorner(idx,1:2);
                bodyLoc(patchInd,:) = bodyXYRightFrontCorner(1,1:2);
            case 4
                location(patchInd,:) = xyRightFrontCorner(idx,3:4);
                bodyLoc(patchInd,:) = bodyXYRightFrontCorner(1,1:2);
            case 5
                location(patchInd,:) = xyRightRearCorner(idx,1:2);
                bodyLoc(patchInd,:) = bodyXYRightRearCorner(1,1:2);
            case 6
                location(patchInd,:) = xyRightRearCorner(idx,3:4);
                bodyLoc(patchInd,:) = bodyXYRightRearCorner(1,1:2);
            case 7
                location(patchInd,:) = xyLeftRearCorner(idx,1:2);
                bodyLoc(patchInd,:) = bodyXYLeftRearCorner(1,1:2);
            case 8
                location(patchInd,:) = xyLeftRearCorner(idx,3:4);
                bodyLoc(patchInd,:) = bodyXYLeftRearCorner(1,1:2);
            case 9
                location(patchInd,:) = xyLeftSide(idx,1:2);
                bodyLoc(patchInd,:) = bodyXYLeftSide(idx,1:2);
            case 10
                location(patchInd,:) = xyLeftSide(idx,3:4);
                bodyLoc(patchInd,:) = bodyXYLeftSide(idx,3:4);
            case 11
                location(patchInd,:) = xyRightSide(idx,1:2);
                bodyLoc(patchInd,:) = bodyXYRightSide(idx,1:2);
            case 12
                location(patchInd,:) = xyRightSide(idx,3:4);
                bodyLoc(patchInd,:) = bodyXYRightSide(idx,3:4);
            case 13
                location(patchInd,:) = xyFront(idx,1:2);
                bodyLoc(patchInd,:) = bodyXYFront(idx,1:2);
            case 14
                location(patchInd,:) = xyFront(idx,3:4);
                bodyLoc(patchInd,:) = bodyXYFront(idx,3:4);
        end
        
        % With the location set, determine the time required to reach the
        % location
        time(patchInd) = angle(patchInd)/v0;
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