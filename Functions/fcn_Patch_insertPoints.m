function patchStruct = fcn_Patch_insertPoints(patchStruct,pointArray)
% fcn_Patch_insertPoints
% Inserts a set of points into a patch structure, preserving the CCCW
% ordering of the points such that a patch graphics object can be plotted
%
% FORMAT:
%
%       updatedPatch = fcn_Patch_insertPoints(patchStruct,pointArray)
%
% INPUTS:
%
%      patchStruct: a structure containing subfields of X and Y coordinates
%      in the following form
%           patchArray{i_patch}.X
%           patchArray{i_patch}.Y
%      Note that i_patch denotes an array of patch structures. Each
%      structure will be plotted separately.
%      pointArray: an N x 2 matrix with column-wise X and Y coordinates
%                  of points to insert into the patch
%
% OUTPUTS:
%
%      patchStruct: the updated patch structure with the points added
%
% DEPENDENCIES:
%
%      fcn_Patch_determineAABB
%
%      ## NOT CURRENTLY USED: fcn_Patch_checkInputsToFunctions
%
% EXAMPLES:
%
%       See the script: script_test_fcn_Patch_insertPoints.m for a full test
%       suite.
%
% This function was written by C. Beal
% Questions or comments? cbeal@bucknell.edu

% Revision history:
%     2022_01_25
%     -- wrote the code

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
    if nargin ~= 2
        error('Incorrect number of input arguments')
    end
    
    % Check to see that there were any points provided
    if isempty(pointArray)
        warning('Empty array of points to add, nothing to do.');
        return
    end
    
    % Check the point array input
    if size(pointArray,2) ~= 2
        error('Wrong format for point array input')
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
NumPoints = size(pointArray,1);


% The first point will be the closest point, the highest index will be the
% one furthest around counter-clockwise, and the lowest index will be the
% most clockwise point. We want to find the most counter-clockwise
% half-space in which the test point fits.

% define the dot product for use in testing half-space membership
dot = @(x,y) x(:)'*y(:);

for i = 1:NumPoints
    if flag_do_debug
        % Address the proper figure
        figure(11);
        clf
        hold on
        hp = patch(patchStruct.pointsX,patchStruct.pointsY,patchStruct.color);
        hpts = plot(patchStruct.pointsX,patchStruct.pointsY,'k*');
        for j = 1:length(patchStruct.pointsX)
            hlab = text(patchStruct.pointsX(j)+1,patchStruct.pointsY(j),num2str(j));
        end
        plot(pointArray(i,1),pointArray(i,2),'rd')
    end
    
    % If there are no existing points, start the points array
    if 0 == length(patchStruct.pointsX)
        patchStruct.pointsX = pointArray(i,1);
        patchStruct.pointsY = pointArray(i,2);
        % If there is only one point, there is no sense of an ordering, so just
        % add the point to the point vectors
    elseif 1 == length(patchStruct.pointsX)
        patchStruct.pointsX = [patchStruct.pointsX(1:end); pointArray(i,1)];
        patchStruct.pointsY = [patchStruct.pointsY(1:end); pointArray(i,2)];
        % If there are two points, need to determine whether the third point
        % is to the right or left of the existing line segment
    elseif 2 == length(patchStruct.pointsX)
        crossProd = cross([patchStruct.pointsX(2,1) - patchStruct.pointsX(1,1);...
            patchStruct.pointsY(2,1) - patchStruct.pointsY(1,1); 0],...
            [pointArray(i,1) - patchStruct.pointsX(1,1);...
            pointArray(i,2) - patchStruct.pointsY(1,1); 0]);
        % If the cross product is negative, the point is to the right of the line segment
        % formed by the existing two points and should go in between the two existing points
        if 1 > crossProd
            patchStruct.pointsX = [patchStruct.pointsX(1); pointArray(i,1); patchStruct.pointsX(2)];
            patchStruct.pointsY = [patchStruct.pointsY(1); pointArray(i,2); patchStruct.pointsY(2)];
            % Otherwise, the point is to the left of the segment and can go after the existing
            % two points
        else
            patchStruct.pointsX = [patchStruct.pointsX; pointArray(i,1)];
            patchStruct.pointsY = [patchStruct.pointsY; pointArray(i,2)];
        end
    else
        % Find the closest patch point to the test point
        [closePt,~] = knnsearch([patchStruct.pointsX patchStruct.pointsY],pointArray(i,:),'nsmethod','exhaustive');
        prevPt = mod(closePt-2,length(patchStruct.pointsX))+1;
        nextPt = mod(closePt,length(patchStruct.pointsX))+1;
        dotProd = dot(pointArray(i,:) - [patchStruct.pointsX(closePt) patchStruct.pointsY(closePt)],[patchStruct.pointsX(nextPt) patchStruct.pointsY(nextPt)] - [patchStruct.pointsX(closePt) patchStruct.pointsY(closePt)]);
        if(0 < dotProd)
            if(1 == nextPt)
                if flag_do_debug
                    fprintf('New point goes after point %d\n',closePt);
                end
                patchStruct.pointsX = [patchStruct.pointsX(1:end); pointArray(i,1)];
                patchStruct.pointsY = [patchStruct.pointsY(1:end); pointArray(i,2)];
            else
                if flag_do_debug
                    fprintf('New point goes between points %d and %d\n',closePt,nextPt);
                end
                patchStruct.pointsX = [patchStruct.pointsX(1:closePt); pointArray(i,1); patchStruct.pointsX(closePt+1:end)];
                patchStruct.pointsY = [patchStruct.pointsY(1:closePt); pointArray(i,2); patchStruct.pointsY(closePt+1:end)];
            end
        else
            if(1 == closePt)
                if flag_do_debug
                    fprintf('New point goes before point %d\n',closePt);
                end
                patchStruct.pointsX = [pointArray(i,1); patchStruct.pointsX(closePt:end)];
                patchStruct.pointsY = [pointArray(i,2); patchStruct.pointsY(closePt:end)];
            else
                if flag_do_debug
                    fprintf('New point goes between points %d and %d\n',prevPt,closePt);
                end
                patchStruct.pointsX = [patchStruct.pointsX(1:prevPt); pointArray(i,1); patchStruct.pointsX(closePt:end)];
                patchStruct.pointsY = [patchStruct.pointsY(1:prevPt); pointArray(i,2); patchStruct.pointsY(closePt:end)];
            end
        end
    end
    if flag_do_debug
        pause;
        % Address the proper figure
        figure(11);
        clf
        hold on
        hp = patch(patchStruct.pointsX,patchStruct.pointsY,patchStruct.color);
        hpts = plot(patchStruct.pointsX,patchStruct.pointsY,'k*');
        for j = 1:length(patchStruct.pointsX)
            hlab = text(patchStruct.pointsX(j)+2,patchStruct.pointsY(j)+1,num2str(j));
        end
    end
end

% Update the axis-aligned bounding box (aabb)
patchStruct = fcn_Patch_determineAABB(patchStruct);


%% Any debugging?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _
%  |  __ \     | |
%  | |  | | ___| |__  _   _  __ _
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_debug
    % Nothing in here yet
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end
end

