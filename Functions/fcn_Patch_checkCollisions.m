function [collFlag,time,angle,location,clearance,bodyLoc] = fcn_Patch_checkCollisions(x0,vehicle,patchArray)
% fcn_Patch_checkCollisions
% Evaluates a circular vehicle trajectory against a series of patches to
% determine whether there will be collisions between the vehicle and the
% outline of any of the patch objects.
%
% ASSUMPTIONS:
%       1) The vehicle moves at constant speed along a circular trajectory.
%       2) The vehicle travels with zero sideslip angle.
%       3) The vehicle is represented by a bounding rectangle.
%       4) The objects are represented by their convex hull.
%       5) A positive trajectory radius indicates CCW travel. Negative
%       indicates CW.
%
% FORMAT:
%
%       [time,angle,location,clearance,bodyLoc] = fcn_Patch_checkCollisions(x0,vehicle,patchArray)
%
% INPUTS:
%
%      x0: a 5 x 1 vector containing the starting (x,y) coordinates of the
%           vehicle, the initial heading, the longitudinal vehicle speed,
%           and the signed trajectory radius in (m,m), radians, m/s, and m.
%      vehicle: a structure containing the vehicle properties, which must
%           include fields a, b, d for the vehicle CG-front axle distance,
%           CG-rear axle distance, and body width, in meters, respectively.
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
%      fcn_Patch_determineAABB
%
%      ## NOT CURRENTLY USED: fcn_Patch_checkInputsToFunctions
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
% 1) look for cases where obstacle is inside vehicle, or goes into vehicle,
% at astart. This would happen, for example, if a barrel gets moved (a
% "ghost barrel") after the zone is mapped and the AV drives through, then
% the data may indicate situations where the AV data shows the AV having
% the object "inside" the vehicle.
%   -- This is easily done outside of this function with polyxpoly()
%2) infinite radii? 
%3) see break case(s)?  
%    - most fixed, but still run into issue with sideswipe before vehicle
%    clears its own body length at the start

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
    
    % Check to see that there were any points provided
    if isempty(patchArray)
        warning('Empty array of objects, nothing to do.');
        return
    end
    
    % Check the trajectory input
    if size(x0,1) ~= 5
        error('Vehicle trajectory missing elements, cannot calculate collisions.')
    end
    
    % Check the vehicle structure input to make sure that the dimensions a,
    % b, and d are all supplied
    if ~all(isfield(vehicle,{'a','b','d'}))
        error('One or more necessary vehicle dimensions missing. Check inputs.')
    end
    if isfield(vehicle,'a') && isempty(vehicle.a)
        error('CG-front axle distance empty.')
    end
    if isfield(vehicle,'b') && isempty(vehicle.a)
        error('CG-rear axle distance empty.')
    end
    if isfield(vehicle,'d') && isempty(vehicle.a)
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
v0 = x0(4);     % Initial (and constant) vehicle speed
R = x0(5);      % Radius of the circular trajectory

% Determine the number of patches to check
Npatches = length(patchArray);

% Calculate center point of vehicle trajectory circle
pc(1) = p0(1) + R*cos(h0+pi/2);
pc(2) = p0(2) + R*sin(h0+pi/2);

% Determine the signed and unsigned bounding radii for all portions of
% the vehicle
Rabs = abs(R);
Rmin = R - vehicle.d/2*sign(R);
RminAbs = abs(Rmin);
Rmax = sign(R)*sqrt((R+vehicle.d/2*sign(R))^2 + max(vehicle.a,vehicle.b)^2);
RmaxAbs = abs(Rmax);

% Find the radii for the front left and front right corners
Rinside = sign(R)*sqrt((R-vehicle.d/2*sign(R))^2 + vehicle.a^2);
RinsideAbs = abs(Rinside);
Routside = sign(R)*sqrt((R+vehicle.d/2*sign(R))^2 + vehicle.a^2);
RoutsideAbs = abs(Routside);

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
    % Determine the point indices that are on the convex hull (CEB: appears
    % not to be necessary)
    % convInds = convhull(patchArray(patchInd).pointsX,patchArray(patchInd).pointsY);
    % Subset the patch points to evaluate just the points on the hull
    xobst = patchArray(patchInd).pointsX;
    yobst = patchArray(patchInd).pointsY;
    % CEB: Determine which of the indices and edges are relevant for this
    % check? If both indices are outside of the annulus *and* the edge does
    % not cross through, it is not relevant here. (Is this correct?)
    % Vertices are relevant if
    % the vertex in the annulus or the next edge crosses through.
    
    % Determine the distance from the center of the vehicle trajectory for
    % all convex hull points of each of the obstacles
    vertexRadii = sqrt((xobst - pc(1)).^2 + (yobst - pc(2)).^2);
    vertexAngles = atan2(yobst - pc(2), xobst - pc(1));
    
    Nobst = size(xobst,1);
    thetaVertex = nan(Nobst,1);
    bodyXYVertex = nan(Nobst,2);
    thetaIFEdge = nan(Nobst,1);
    xyIFEdge = nan(Nobst,2);
    bodyXYIFEdge = nan(Nobst,2);
    thetaOFEdge = nan(Nobst,1);
    xyOFEdge = nan(Nobst,2);
    bodyXYOFEdge = nan(Nobst,2);
    thetaOREdge = nan(Nobst,1);
    xyOREdge = nan(Nobst,2);
    bodyXYOREdge = nan(Nobst,2);
    
    for vertexInd = 1:Nobst
        % First, determine if the vertex will collide with the front of the car
        if vertexRadii(vertexInd) <= RoutsideAbs && vertexRadii(vertexInd) >= RinsideAbs
            % Calculate based on the front of the car
            thetaVertex(vertexInd) = vertexAngles(vertexInd) - sign(R)*asin(vehicle.a/vertexRadii(vertexInd));
            thetaVertex(vertexInd) = rerangeAngles(thetaVertex(vertexInd));
            bodyXYVertex(vertexInd,:) = [vehicle.a R - vertexRadii(vertexInd)];
        elseif vertexRadii(vertexInd) <= RmaxAbs && vertexRadii(vertexInd) >= RminAbs
            % Calculate based on the sides of the car
            if vertexRadii(vertexInd) <= RinsideAbs && vertexRadii(vertexInd) >= RminAbs
                % Calculate based on the inside of the car. The nearest point
                % is alongside the CG, so the collision will always happen
                % ahead of the CG.
                alphaa = sqrt(vertexRadii(vertexInd)^2 - (Rabs-vehicle.d/2)^2);
                % Adjust the angle by the computed distance ahead of the CG.
                thetaVertex(vertexInd) = vertexAngles(vertexInd) - sign(R)*atan(alphaa/(Rabs-vehicle.d/2));
                thetaVertex(vertexInd) = rerangeAngles(thetaVertex(vertexInd));
                bodyXYVertex(vertexInd,:) = [alphaa vehicle.d/2];
            else
                % Calculate based on the outside of the car. If the object
                % cleared the front, it will only hit behind the CG.
                alphab = sqrt(vertexRadii(vertexInd)^2 - (Rabs+vehicle.d/2)^2);
                % Adjust the angle by the computed distance behind the CG.
                thetaVertex(vertexInd) = vertexAngles(vertexInd) + sign(R)*atan(alphab/(Rabs+vehicle.d/2));
                thetaVertex(vertexInd) = rerangeAngles(thetaVertex(vertexInd));
                bodyXYVertex(vertexInd,:) = [-alphab vehicle.d/2];
            end
        end
        % Determine the index of the next vertex to define edges
        nextVertex = mod(vertexInd,Nobst)+1;
        
        % Next, check all of the obstacle edges for "spanning" edges where
        % one or both vertices of the convex hull are outside of the
        % trajectory but the object actually lies across the trajectory.
        
        % Check for intersections with the front of the vehicle
        pa = [xobst(vertexInd) yobst(vertexInd)];
        pb = [xobst(nextVertex) yobst(nextVertex)];
        
        % Compute the alpha associated with the min radius
        alpha = -((pa(1)-pc(1))*(pb(1)-pa(1)) + (pa(2)-pc(2))*(pb(2)-pa(2)))/((pb(1)-pa(1))^2 + (pb(2)-pa(2))^2);
        % Bound alpha between 0 and 1;
        % alpha = max(0,min(alpha,1));
        
        % Check to see if there is a segment that spans the outer radius
        if(vertexRadii(vertexInd) > RmaxAbs && vertexRadii(nextVertex) > RmaxAbs && alpha > 0 && alpha < 1)
            % Handle this case, somehow. Best guess is to somehow divide
            % the segment into two and use the existing code. However, this
            % requires either dividing the segment before this inside loop
            % or shifting the vertices and getting the 'for' loop to expand
            % the number of iterations accordingly
            pd(1) = pa(1) + alpha*(pb(1)-pa(1));
            pd(2) = pa(2) + alpha*(pb(2)-pa(2));
            RsegmentMin = sqrt((pd(1)-pc(1))^2 + (pd(2)-pc(2))^2);
            % If the resulting minimum radius is less than the outer radius
            % of the vehicle path, the edge spans (since the end points are
            % outside based on the outer conditional)
            if RsegmentMin < RmaxAbs
                % Split the edge into two parts and test both of them for
                % intersections with the outer circle
                [thetaOFEdge(vertexInd),xyOFEdge(vertexInd,:)] = intersectEdgeWithCircle(pa,pd,pc,RoutsideAbs);
                [thetaOFTest,xyOFTest] = intersectEdgeWithCircle(pb,pd,pc,RoutsideAbs);
                if thetaOFTest < thetaOFEdge(vertexInd)
                    thetaOFEdge(vertexInd) = thetaOFTest;
                    xyOFEdge(vertexInd,:) = xyOFTest;
                end
                thetaOFEdge(vertexInd) = thetaOFEdge(vertexInd) - sign(R)*asin(vehicle.a/RoutsideAbs);
                bodyXYOFEdge(vertexInd,:) = [vehicle.a vehicle.d/2];
                % Split the edge into two parts and test both of them for
                % intersections with the outer circle
                [thetaOREdge(vertexInd),xyOREdge(vertexInd,:)] = intersectEdgeWithCircle(pa,pd,pc,RmaxAbs);
                [thetaORTest,xyORTest] = intersectEdgeWithCircle(pb,pd,pc,RmaxAbs);
                if thetaORTest < thetaOREdge(vertexInd)
                    thetaOREdge(vertexInd) = thetaORTest;
                    xyOREdge(vertexInd,:) = xyORTest;
                end
                thetaOREdge(vertexInd) = thetaOREdge(vertexInd) + sign(R)*asin(vehicle.b/RmaxAbs);
                bodyXYOREdge(vertexInd,:) = [-vehicle.b vehicle.d/2];
            end
        else
            [thetaOFEdge(vertexInd),xyOFEdge(vertexInd,:)] = intersectEdgeWithCircle(pa,pb,pc,RoutsideAbs);
            thetaOFEdge(vertexInd) = thetaOFEdge(vertexInd) - sign(R)*asin(vehicle.a/RoutsideAbs);
            bodyXYOFEdge(vertexInd,:) = [vehicle.a vehicle.d/2];
            [thetaOREdge(vertexInd),xyOREdge(vertexInd,:)] = intersectEdgeWithCircle(pa,pb,pc,RmaxAbs);
            thetaOREdge(vertexInd) = thetaOREdge(vertexInd) + sign(R)*asin(vehicle.b/RmaxAbs);
            bodyXYOREdge(vertexInd,:) = [-vehicle.b vehicle.d/2];
        end
        [thetaIFEdge(vertexInd),xyIFEdge(vertexInd,:)] = intersectEdgeWithCircle(pa,pb,pc,RinsideAbs);
        thetaIFEdge(vertexInd) = thetaIFEdge(vertexInd) - sign(R)*asin(vehicle.a/RinsideAbs);
        bodyXYIFEdge(vertexInd,:) = [vehicle.a -vehicle.d/2];
        
    end
    % Check for a case where no intersections were calculated
    if all([isnan(thetaVertex); isnan(thetaIFEdge); isnan(thetaOFEdge); isnan(thetaOREdge)])
        collFlag(patchInd) = 0;
        
        % Determine the nearest two vertices
        [nearestOuterClearance,nearestOuters] = mink(min(inf(Nobst,1),vertexRadii - RmaxAbs),2);
        [nearestInnerClearance,nearestInners] = mink(min(inf(Nobst,1),RminAbs - vertexRadii),2);
        if max(0,nearestOuterClearance) > max(nearestInnerClearance)
            pa = [xobst(nearestOuters(1)) yobst(nearestOuters(1))];
            pb = [xobst(nearestOuters(2)) yobst(nearestOuters(2))];
            theta_offset = sign(R)*asin(vehicle.b/Rabs);
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
            time(patchInd) = rerangeAngles(angle(patchInd) - h0 + pi/2)*Rabs/v0;
        else
            time(patchInd) = rerangeAngles(h0 + pi/2 - angle(patchInd))*Rabs/v0;
        end
        % No collision, so no body collision location
        bodyLoc(patchInd,:) = [NaN NaN];
        
        % Now, find the minimum angular location for the patch object in the
        % vehicle travel direction (indicated by the sign of R)
    elseif R >= 0
        collFlag(patchInd) = 1;
        travelOffset = -h0 + pi/2;
        [minVertex,minVertexInd] = nanmin(rerangeAngles(thetaVertex+travelOffset));
        [minIFEdge,minIFEdgeInd] = nanmin(rerangeAngles(thetaIFEdge+travelOffset));
        [minOFEdge,minOFEdgeInd] = nanmin(rerangeAngles(thetaOFEdge+travelOffset));
        [minOREdge,minOREdgeInd] = nanmin(rerangeAngles(thetaOREdge+travelOffset));
        if isnan(minVertex)
            minVertex = inf;
        end
        if isnan(minIFEdge)
            minIFEdge = inf;
        end
        if isnan(minOFEdge)
            minOFEdge = inf;
        end
        if isnan(minOREdge)
            minOREdge = inf;
        end
        if minVertex < minIFEdge && minVertex < minOFEdge && minVertex < minOREdge
            angle(patchInd) = thetaVertex(minVertexInd);
            location(patchInd,:) = [xobst(minVertexInd) yobst(minVertexInd)];
            bodyLoc(patchInd,:) = bodyXYVertex(minVertexInd,:);
        elseif minIFEdge < minVertex && minIFEdge < minOFEdge && minIFEdge < minOREdge
            angle(patchInd) = thetaIFEdge(minIFEdgeInd);
            location(patchInd,:) = [xyIFEdge(minIFEdgeInd,1) xyIFEdge(minIFEdgeInd,2)];
            bodyLoc(patchInd,:) = bodyXYIFEdge(minIFEdgeInd,:);
        elseif minOFEdge < minVertex && minOFEdge < minIFEdge && minOFEdge < minOREdge
            angle(patchInd) = thetaOFEdge(minOFEdgeInd);
            location(patchInd,:) = [xyOFEdge(minOFEdgeInd,1) xyOFEdge(minOFEdgeInd,2)];
            bodyLoc(patchInd,:) = bodyXYOFEdge(minOFEdgeInd,:);
        else
            angle(patchInd) = thetaOREdge(minOREdgeInd);
            location(patchInd,:) = [xyOREdge(minOREdgeInd,1) xyOREdge(minOREdgeInd,2)];
            bodyLoc(patchInd,:) = bodyXYOREdge(minOREdgeInd,:);
        end
        % With the location set, determine the time required to reach the
        % location
        time(patchInd) = rerangeAngles(angle(patchInd) - h0 + pi/2)*Rabs/v0;
    else
        collFlag(patchInd) = 1;
        % Shift the angles into the negative range for comparison with the
        % max function (this isn't correct)
        travelOffset = h0 + pi/2;
        [minVertex,minVertexInd] = nanmin(rerangeAngles(travelOffset - thetaVertex));
        [minIFEdge,minIFEdgeInd] = nanmin(rerangeAngles(travelOffset - thetaIFEdge));
        [minOFEdge,minOFEdgeInd] = nanmin(rerangeAngles(travelOffset - thetaOFEdge));
        [minOREdge,minOREdgeInd] = nanmin(rerangeAngles(travelOffset - thetaOREdge));
        if isnan(minVertex)
            minVertex = inf;
        end
        if isnan(minIFEdge)
            minIFEdge = inf;
        end
        if isnan(minOFEdge)
            minOFEdge = inf;
        end
        if isnan(minOREdge)
            minOREdge = inf;
        end
        if minVertex < minIFEdge && minVertex < minOFEdge && minVertex < minOREdge
            angle(patchInd) = thetaVertex(minVertexInd);
            location(patchInd,:) = [xobst(minVertexInd) yobst(minVertexInd)];
            bodyLoc(patchInd,:) = bodyXYVertex(minVertexInd,:);
        elseif minIFEdge < minVertex && minIFEdge < minOFEdge && minIFEdge < minOREdge
            angle(patchInd) = thetaIFEdge(minIFEdgeInd);
            location(patchInd,:) = [xyIFEdge(minIFEdgeInd,1) xyIFEdge(minIFEdgeInd,2)];
            bodyLoc(patchInd,:) = bodyXYIFEdge(minIFEdgeInd,:);
        elseif minOFEdge < minVertex && minOFEdge < minIFEdge && minOFEdge < minOREdge
            angle(patchInd) = thetaOFEdge(minOFEdgeInd);
            location(patchInd,:) = [xyOFEdge(minOFEdgeInd,1) xyOFEdge(minOFEdgeInd,2)];
            bodyLoc(patchInd,:) = bodyXYOFEdge(minOFEdgeInd,:);
        else
            angle(patchInd) = thetaOREdge(minOREdgeInd);
            location(patchInd,:) = [xyOREdge(minOREdgeInd,1) xyOREdge(minOREdgeInd,2)];
            bodyLoc(patchInd,:) = bodyXYOREdge(minOREdgeInd,:);
        end
        % With the location set, determine the time required to reach the
        % location
        time(patchInd) = rerangeAngles(h0 + pi/2 - angle(patchInd))*Rabs/v0;
        clearance(patchInd) = NaN;
    end
end

end

% Function to re-range angles between 0 and 2*pi
function outAngles = rerangeAngles(inAngles)
outAngles = inAngles;
while any(outAngles > 2*pi)
    outAngles(outAngles > 2*pi) = outAngles(outAngles > 2*pi) - 2*pi;
end
while any(outAngles < 0)
    outAngles(outAngles < 0) = outAngles(outAngles < 0) + 2*pi;
end
end

% Function to intersect an edge defined as the segment between two points
% pa->[x y] and pb->[x y] with a circle of radius R and centered at point
% pc->[x y]. The return values are the angle of intersection and the point
% pint->[x y] where the intersection occurs. This does not handle cases
% where the edge could be larger than the circle itself and thus have two
% intersections
function [intAngle,intPoint] = intersectEdgeWithCircle(pa,pb,pc,R)

% First check to see if an intersection is possible
Rabs = abs(R);
if (norm(pa-pc) < Rabs && norm(pb-pc) < Rabs) || (norm(pa-pc) > Rabs && norm(pb-pc) > Rabs)
    % If no intersection (or two intersections) possible, return NaNs
    intAngle = NaN;
    intPoint = [NaN NaN];
else
    % Calculate the coefficients of the quadratic equation for alpha
    % (the scalar distance along the line segment)
    quadCoefs(1) = (pb(1) - pa(1))^2 + (pb(2) - pa(2))^2;
    quadCoefs(2) = 2*((pa(1)-pc(1))*(pb(1) - pa(1)) + (pa(2)-pc(2))*(pb(2) - pa(2)));
    quadCoefs(3) = (pa(1)-pc(1))^2 + (pa(2)-pc(2))^2 - R^2;
    % Check the conditions on the coefficients to obtain alpha between
    % 0 and 1 and solve for alpha. The quadratic coefficients 1 and 3
    % should always be positive. The discriminant should always be
    % positive or zero as well, since we eliminate cases where there is
    % no intersection or multiple intersections. Thus, we only need to
    % worry about whether to use the root from the + or - operation in
    % the quadratic formula. We can determine this by testing to see if
    % the first term in the quadratic formula is greater than 1. If so,
    % we need to use the root with the negative sign. Otherwise, we use
    % the root with the positive sign.
    %     if -quadCoefs(2)/(2*quadCoefs(1)) < 0
    %         alpha = (-quadCoefs(2) - sqrt(quadCoefs(2)^2 - 4*quadCoefs(1)*quadCoefs(3)))/(2*quadCoefs(1));
    %     else
    %         alpha = (-quadCoefs(2) + sqrt(quadCoefs(2)^2 - 4*quadCoefs(1)*quadCoefs(3)))/(2*quadCoefs(1));
    %     end
    % Calculate the additive root for alpha
    alpha = (-quadCoefs(2) + sqrt(quadCoefs(2)^2 - 4*quadCoefs(1)*quadCoefs(3)))/(2*quadCoefs(1));
    % Check to see if the solution should come from the subtractive root
    % and adjust if necessary
    if alpha > 1
        alpha = (-quadCoefs(2) - sqrt(quadCoefs(2)^2 - 4*quadCoefs(1)*quadCoefs(3)))/(2*quadCoefs(1));
    end
    % Confirm that we got a good solution. (This shouldn't fail with
    % the previous logic, but it's here as extra validation.)
    if alpha > 1 || alpha < 0
        error('Calculation error, out of range');
    end
    % Calculate the associated circular angle of intersection from the
    % scalar distance along the line segment (alpha)
    intAngle = atan2(pa(2) + alpha*(pb(2)-pa(2)) - pc(2),...
        pa(1) + alpha*(pb(1)-pa(1)) - pc(1));
    % Normalize the range of the intersection angle to [0,2*pi]
    intAngle = rerangeAngles(intAngle);
    % Compute the actual point of intersection
    intPoint(1) = pa(1) + alpha*(pb(1)-pa(1));
    intPoint(2) = pa(2) + alpha*(pb(2)-pa(2));
end
end
