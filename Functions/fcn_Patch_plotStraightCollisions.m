function fcn_Patch_plotStraightCollisions(collFlags,collLoc,clearance,bodyCollLoc,vehicle,h0,varargin)
% fcn_Patch_plotPatch
% Plots a visual representation of the objects in a patch structure array
%
% FORMAT: 
%
%       fcn_Patch_plotStraightCollisions(collFlags,collLoc,clearance,bodyCollLoc,vehicle,h0,(fig_num))
%
% INPUTS:
%
%      collFlags: an N x 1 vector of flags denoting whether there is a
%           collision with each of the objects
%      collLoc: a N x 2 vector of collision locations, where N is the
%           number of patch objects and the columns are the x and y
%           coordinates of the collision. Elements of the location matrix
%           will be set to NaN if there is no overlap with the vehicle
%           path.
%       clearance: an N x 1 vector of minimum clearance distances between
%           the vehicle and the patch objects, where N is the number of
%           patch objects. Elements of the clearance vector will be set to
%           NaN if there is a collision with the object.
%       bodyCollLoc: an N x 2 vector of collision locations in vehicle body
%           fixed coordinates, where N is the number of patch objects and
%           the columns are the x and y coordinates of the collision.
%           Elements of the location matrix will be set to NaN if there is
%           no overlap with the vehicle path
%       vehicle: [1x1] cell struct with vehicle properties
%       h0: [1x1] scalar, representing initial heading of the vehicle
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
%       See the script: script_test_fcn_Patch_plotStraightCollisions.m for a full test
%       suite. 
%
% This function was written by Shashank Garikipati

% Revision history:
%     2022_07_06
%     -- wrote the code

flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking
flag_mediumRes = 1;

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
if 6 < nargin
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

    plot(bodyLoc([1:end 1],1),bodyLoc([1:end 1],2),'b-','linewidth',1);

    if 1 == flag_mediumRes
    
        body1Loc(1,1) = collLoc(collInd,1) + (vehicle.df - plotOffset(1))*cos(h0) - (((vehicle.w/2-(1/2)*vehicle.w/2))-plotOffset(2))*sin(h0);
        body1Loc(1,2) = collLoc(collInd,2) + (vehicle.df - plotOffset(1))*sin(h0) + (((vehicle.w/2-(1/2)*vehicle.w/2))-plotOffset(2))*cos(h0);
        body1Loc(2,1) = collLoc(collInd,1) + ((vehicle.df-(1/3.5)*vehicle.df) - plotOffset(1))*cos(h0) - (vehicle.w/2-plotOffset(2))*sin(h0);
        body1Loc(2,2) = collLoc(collInd,2) + ((vehicle.df-(1/3.5)*vehicle.df) - plotOffset(1))*sin(h0) + (vehicle.w/2-plotOffset(2))*cos(h0);
        body1Loc(3,1) = collLoc(collInd,1) + ((-vehicle.dr/2-(1/8)*vehicle.dr/2) - plotOffset(1))*cos(h0) - (vehicle.w/2-plotOffset(2))*sin(h0);
        body1Loc(3,2) = collLoc(collInd,2) + ((-vehicle.dr/2-(1/8)*vehicle.dr/2) - plotOffset(1))*sin(h0) + (vehicle.w/2-plotOffset(2))*cos(h0);
        body1Loc(4,1) = collLoc(collInd,1) + (-vehicle.dr - plotOffset(1))*cos(h0) - ((vehicle.w/2-(1/2)*vehicle.w/2)-plotOffset(2))*sin(h0);
        body1Loc(4,2) = collLoc(collInd,2) + (-vehicle.dr - plotOffset(1))*sin(h0) + ((vehicle.w/2-(1/2)*vehicle.w/2)-plotOffset(2))*cos(h0);
        
        body1Loc(5,1) = collLoc(collInd,1) + (-vehicle.dr - plotOffset(1))*cos(h0) - (((-vehicle.w/2+(1/2)*vehicle.w/2))-plotOffset(2))*sin(h0);
        body1Loc(5,2) = collLoc(collInd,2) + (-vehicle.dr - plotOffset(1))*sin(h0) + (((-vehicle.w/2+(1/2)*vehicle.w/2))-plotOffset(2))*cos(h0);
        body1Loc(6,1) = collLoc(collInd,1) + ((-vehicle.dr/2-(1/8)*vehicle.dr/2) - plotOffset(1))*cos(h0) - (-vehicle.w/2-plotOffset(2))*sin(h0);
        body1Loc(6,2) = collLoc(collInd,2) + ((-vehicle.dr/2-(1/8)*vehicle.dr/2) - plotOffset(1))*sin(h0) + (-vehicle.w/2-plotOffset(2))*cos(h0);
        body1Loc(7,1) = collLoc(collInd,1) + ((vehicle.df-(1/3.5)*vehicle.df) - plotOffset(1))*cos(h0) - (-vehicle.w/2-plotOffset(2))*sin(h0);
        body1Loc(7,2) = collLoc(collInd,2) + ((vehicle.df-(1/3.5)*vehicle.df) - plotOffset(1))*sin(h0) + (-vehicle.w/2-plotOffset(2))*cos(h0);
        body1Loc(8,1) = collLoc(collInd,1) + (vehicle.df - plotOffset(1))*cos(h0) - (((-vehicle.w/2+(1/2)*vehicle.w/2))-plotOffset(2))*sin(h0);
        body1Loc(8,2) = collLoc(collInd,2) + (vehicle.df - plotOffset(1))*sin(h0) + (((-vehicle.w/2+(1/2)*vehicle.w/2))-plotOffset(2))*cos(h0);
        plot(body1Loc([1:end 1],1),body1Loc([1:end 1],2),'k','linewidth',1);
    
        body2Loc(1,1) = collLoc(collInd,1) + ((vehicle.df-(1/3)*vehicle.df) - plotOffset(1))*cos(h0) - (vehicle.w/2-plotOffset(2))*sin(h0);
        body2Loc(1,2) = collLoc(collInd,2) + ((vehicle.df-(1/3)*vehicle.df) - plotOffset(1))*sin(h0) + (vehicle.w/2-plotOffset(2))*cos(h0);
        body2Loc(2,1) = collLoc(collInd,1) + ((vehicle.df-(1/3.8)*vehicle.df) - plotOffset(1))*cos(h0) - ((vehicle.w/2-(1/2)*vehicle.w/2)-plotOffset(2))*sin(h0);
        body2Loc(2,2) = collLoc(collInd,2) + ((vehicle.df-(1/3.8)*vehicle.df) - plotOffset(1))*sin(h0) + ((vehicle.w/2-(1/2)*vehicle.w/2)-plotOffset(2))*cos(h0);
        body2Loc(3,1) = collLoc(collInd,1) + ((vehicle.df-(1/3.8)*vehicle.df) - plotOffset(1))*cos(h0) - ((-vehicle.w/2+(1/2)*vehicle.w/2)-plotOffset(2))*sin(h0);
        body2Loc(3,2) = collLoc(collInd,2) + ((vehicle.df-(1/3.8)*vehicle.df) - plotOffset(1))*sin(h0) + ((-vehicle.w/2+(1/2)*vehicle.w/2)-plotOffset(2))*cos(h0);
        body2Loc(4,1) = collLoc(collInd,1) + ((vehicle.df-(1/3)*vehicle.df) - plotOffset(1))*cos(h0) - ((-vehicle.w/2)-plotOffset(2))*sin(h0);
        body2Loc(4,2) = collLoc(collInd,2) + ((vehicle.df-(1/3)*vehicle.df) - plotOffset(1))*sin(h0) + ((-vehicle.w/2)-plotOffset(2))*cos(h0);
    
        body2Loc(5,1) = collLoc(collInd,1) + ((vehicle.df-(1/1.8)*vehicle.df) - plotOffset(1))*cos(h0) - ((-vehicle.w/2+(1/4)*vehicle.w/2)-plotOffset(2))*sin(h0);
        body2Loc(5,2) = collLoc(collInd,2) + ((vehicle.df-(1/1.8)*vehicle.df) - plotOffset(1))*sin(h0) + ((-vehicle.w/2+(1/4)*vehicle.w/2)-plotOffset(2))*cos(h0);
        body2Loc(6,1) = collLoc(collInd,1) + ((vehicle.df-(1/1.8)*vehicle.df) - plotOffset(1))*cos(h0) - ((vehicle.w/2-(1/4)*vehicle.w/2)-plotOffset(2))*sin(h0);
        body2Loc(6,2) = collLoc(collInd,2) + ((vehicle.df-(1/1.8)*vehicle.df) - plotOffset(1))*sin(h0) + ((vehicle.w/2-(1/4)*vehicle.w/2)-plotOffset(2))*cos(h0);
    
        plot(body2Loc([1:end 1],1),body2Loc([1:end 1],2),'k','linewidth',1);
    
        body3Loc(1,1) = collLoc(collInd,1) + ((-vehicle.dr/6) - plotOffset(1))*cos(h0) - ((vehicle.w/2-(1.5/4)*vehicle.w/2)-plotOffset(2))*sin(h0);
        body3Loc(1,2) = collLoc(collInd,2) + ((-vehicle.dr/6) - plotOffset(1))*sin(h0) + ((vehicle.w/2-(1.5/4)*vehicle.w/2)-plotOffset(2))*cos(h0);
        body3Loc(2,1) = collLoc(collInd,1) + ((-vehicle.dr/6) - plotOffset(1))*cos(h0) - (-(vehicle.w/2-(1.5/4)*vehicle.w/2)-plotOffset(2))*sin(h0);
        body3Loc(2,2) = collLoc(collInd,2) + ((-vehicle.dr/6) - plotOffset(1))*sin(h0) + (-(vehicle.w/2-(1.5/4)*vehicle.w/2)-plotOffset(2))*cos(h0);
        body3Loc(3,1) = collLoc(collInd,1) + ((-vehicle.dr/1.5) - plotOffset(1))*cos(h0) - (-(vehicle.w/2-(1.5/4)*vehicle.w/2)-plotOffset(2))*sin(h0);
        body3Loc(3,2) = collLoc(collInd,2) + ((-vehicle.dr/1.5) - plotOffset(1))*sin(h0) + (-(vehicle.w/2-(1.5/4)*vehicle.w/2)-plotOffset(2))*cos(h0);
        body3Loc(4,1) = collLoc(collInd,1) + ((-vehicle.dr/1.4) - plotOffset(1))*cos(h0) - (-plotOffset(2))*sin(h0);
        body3Loc(4,2) = collLoc(collInd,2) + ((-vehicle.dr/1.4) - plotOffset(1))*sin(h0) + (-plotOffset(2))*cos(h0);
        body3Loc(5,1) = collLoc(collInd,1) + ((-vehicle.dr/1.5) - plotOffset(1))*cos(h0) - ((vehicle.w/2-(1.5/4)*vehicle.w/2)-plotOffset(2))*sin(h0);
        body3Loc(5,2) = collLoc(collInd,2) + ((-vehicle.dr/1.5) - plotOffset(1))*sin(h0) + ((vehicle.w/2-(1.5/4)*vehicle.w/2)-plotOffset(2))*cos(h0);

        plot(body3Loc([1:end 1],1),body3Loc([1:end 1],2),'k','linewidth',1);

        if 0 == collFlags(:,collInd)
            a = sqrt((body1Loc(1,1) - collLoc(collInd,1))^2 + (body1Loc(1,2) - collLoc(collInd,2))^2);
            b = sqrt((body1Loc(4,1) - collLoc(collInd,1))^2 + (body1Loc(4,2) - collLoc(collInd,2))^2);
            if b > a
                quiver(body1Loc(1,1),body1Loc(1,2) ,collLoc(collInd,1) - body1Loc(1,1),collLoc(collInd,2) - body1Loc(1,2),0,'k','linewidth',1);
            else
                quiver(body1Loc(4,1),body1Loc(4,2) ,collLoc(collInd,1) - body1Loc(4,1),collLoc(collInd,2) - body1Loc(4,2),0,'k','linewidth',1);
            end
        end
    end

    % Plot the collision or closest clearance point
    plot(collLoc(collInd,1),collLoc(collInd,2),'r*')
end



