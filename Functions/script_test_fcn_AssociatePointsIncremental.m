% Script to process data sets and test for congruency

clear all
close all

% define the dot product
dot = @(x,y) x(:)'*y(:)

% Generate three data sets via interactive plot
% figure(9);
% hold on
% axis equal
% Set the title header
% UserData.title_header = sprintf('Dataset %.0d of %.0d',1,3);
% % Save the results
% set(gcf,'UserData',UserData);
% inputData = fcn_Dataset_fillSetViaUserInputs(9);
% xyData{1} = inputData;
% % Set the title header
% UserData.title_header = sprintf('Dataset %.0d of %.0d',2,3);
% % Save the results
% set(gcf,'UserData',UserData);
% inputData = fcn_Dataset_fillSetViaUserInputs(9);
% xyData{2} = inputData;
% % Set the title header
% UserData.title_header = sprintf('Dataset %.0d of %.0d',3,3);
% % Save the results
% set(gcf,'UserData',UserData);
% inputData = fcn_Dataset_fillSetViaUserInputs(9);
% xyData{3} = inputData;

% Assign a region of interest around the rover
ROIradius = 20;

% Assign a radial standard deviation value
sigma = 3; % meters

% % Load up some data (simple xy points for now)
% %xyData = fcn_Dataset_fillSampleSets;
load testDatasetVehicle.mat

%% Plot the data (known to be two sets for now)
figure(1)
clf
hold on
for i = 1:length(xyData{1})
    plot(xyData{1}(i,1),xyData{1}(i,2),'b.','Markersize',20,'Markerfacecolor','white');
    %text(xyData{1}(i,1)+2,xyData{1}(i,2),sprintf('%d',i))
end
for i = 1:length(xyData{2})
    plot(xyData{2}(i,1),xyData{2}(i,2),'r.','Markersize',20,'Markerfacecolor','white');
    %text(xyData{2}(i,1)+2,xyData{2}(i,2),sprintf('%d',i))
end
for i = 1:length(xyData{3})
    plot(xyData{3}(i,1),xyData{3}(i,2),'kd','Markersize',10,'Markerfacecolor','white');
end

% Define a red-yellow-green color map for the quality of the data
cMap = interp1([0;1],[0 1 0; 1 0 0],linspace(0,1,256));


%% Go ROI by ROI through the map (as if a vehicle was driving through the map)
for locInd = 1:length(xyData{3})
    % Set up for another round of drawing on the figure
    figure(8)
    clf
    hold on
    axis equal
    
    % Define a vehicle travel direction vector
    if locInd > 1
        travelVec = xyData{3}(locInd,:) - xyData{3}(locInd-1,:);
    else
        travelVec = xyData{3}(2,:) - xyData{3}(1,:);
    end
    
    % Determine the map points in a circular ROI around the vehicle
    mapIdxROI = rangesearch(xyData{1},xyData{3}(locInd,:),ROIradius);
    % Define counting variables to loop through the ROI indices
    counter = 1;
    totalElems = length(mapIdxROI{:});
    % Loop through the map ROI indices and check if they are in front of
    % the vehicle. If not, discard them and shorten the vector of indices.
    while counter <= totalElems
        % Check the dot product of the travel vector and the vector to the
        % data point from the vehicle
        if 0 > dot(travelVec,xyData{1}(mapIdxROI{:}(counter),:) - xyData{3}(locInd,:))
            % If the dot product is negative, the point is behind the
            % vehicle and can be removed
            mapIdxROI = {[mapIdxROI{:}(1:counter-1) mapIdxROI{:}(counter+1:end)]};
            % compensate for having removed one entry
            totalElems = totalElems - 1;
        else
            % Only increment the counter if no points were removed. If
            % points were removed, the next one will have been shifted into
            % the current vector location indexed by the counter and there
            % is no need to increment.
            counter = counter + 1;
        end
    end
    % Determine the rover points in a circular ROI around the vehicle
    roverIdxROI = rangesearch(xyData{2},xyData{3}(locInd,:),ROIradius);
    % Define counting variables to loop through the ROI indices
    counter = 1;
    totalElems = length(roverIdxROI{:});
    % Loop through the rover ROI indices and check if they are in front of
    % the vehicle. If not, discard them and shorten the vector of indices.
    while counter <= totalElems
        % Check the dot product of the travel vector and the vector to the
        % data point from the vehicle
        if 0 > dot(travelVec,xyData{2}(roverIdxROI{:}(counter),:) - xyData{3}(locInd,:))
            % If the dot product is negative, the point is behind the
            % vehicle and can be removed
            roverIdxROI = {[roverIdxROI{:}(1:counter-1) roverIdxROI{:}(counter+1:end)]};
            % compensate for having removed one entry
            totalElems = totalElems - 1;
        else
            % Only increment the counter if no points were removed. If
            % points were removed, the next one will have been shifted into
            % the current vector location indexed by the counter and there
            % is no need to increment.
            counter = counter + 1;
        end
    end
    
    % Subset the map and rover data to handle only the relevant
    % semicircular ROI
    mapDataROI = xyData{1}(mapIdxROI{:},:);
    roverDataROI = xyData{2}(roverIdxROI{:},:);
    
    % Binary compatibility analysis (mutually nearest neighbor points)
    % First, associate data points in the sets using a nearest neighbor search
    [idx1to2,dist1] = knnsearch(roverDataROI,mapDataROI);
    [idx2to1,dist2] = knnsearch(mapDataROI,roverDataROI);
    
    Nmatch = length(idx1to2);
    
    if Nmatch > 0
        % Find mutual matches between nearest neighbors. These points are "most" in
        % agreement in that they are mutual nearest neighbor matches
        mutual = nan(Nmatch,1);
        % Run through the indices of matches from rover data to the map
        for i = 1:length(idx1to2)
            % If there is a mutual match
            if(idx2to1(idx1to2(i)) == i)
                % Note the match in a result vector
                mutual(i) = 1;
            else
                mutual(i) = 0;
            end
        end
        
        % Determine the binary compatibility of the data sets (percentage of points
        % that are mutually nearest neighbors)
        binMetric = sum(mutual)/length(mutual);
        fprintf("%.2f percent of map features matched in rover data set\n",100*binMetric);
        % Determine the rms error between mutually compatible points
        rmsError = sqrt(mean(dist1(mutual == 1).^2));
        fprintf("RMS error for matched features is %.2f units\n",rmsError);
    else
        fprintf("No matched features in the ROI\n");
        mutual = [];
    end
    
    % Determine the vehicle heading angle from the travel vector
    startAngle = atan2(travelVec(2),travelVec(1));
    % Create a semicircle centered around the vehicle heading angle
    semicrc = ROIradius.*[cos((-pi/2:0.1:pi/2) + startAngle); sin((-pi/2:0.1:pi/2) + startAngle)];
    if Nmatch > 0
        patch(semicrc(1,:) + xyData{3}(locInd,1), semicrc(2,:) + xyData{3}(locInd,2), interp1(linspace(0.5,2,256)',cMap,rmsError,'nearest','extrap'))
    else
        patch(semicrc(1,:) + xyData{3}(locInd,1), semicrc(2,:) + xyData{3}(locInd,2), [0.4 0.4 0.4])    
    end
    for i = 1:length(xyData{1})
        plot(xyData{1}(i,1),xyData{1}(i,2),'bo','Markersize',10,'Markerfacecolor','white');
        %text(xyData{1}(i,1)+2,xyData{1}(i,2),sprintf('%d',i))
    end
    for i = 1:length(xyData{2})
        plot(xyData{2}(i,1),xyData{2}(i,2),'ro','Markersize',10,'Markerfacecolor','white');
        %text(xyData{2}(i,1)+2,xyData{2}(i,2),sprintf('%d',i))
    end
    plot(xyData{3}(locInd,1),xyData{3}(locInd,2),'kd','Markersize',10,'MarkerFaceColor','black');
    for i = 1:size(mapDataROI,1)
        plot(mapDataROI(i,1),mapDataROI(i,2),'bo','Markersize',10,'Markerfacecolor','b');
        %text(xyData{1}(i,1)+2,xyData{1}(i,2),sprintf('%d',i))
    end
    for i = 1:size(roverDataROI,1)
        plot(roverDataROI(i,1),roverDataROI(i,2),'ro','Markersize',10,'Markerfacecolor','r');
        %text(xyData{2}(i,1)+2,xyData{2}(i,2),sprintf('%d',i))
    end
    plot(xyData{3}(locInd,1)+[0 travelVec(1)],xyData{3}(locInd,2)+[0 travelVec(2)],'k-.')
    % Indicate on the plot that these points match
    for i = 1:length(mutual)
        if 1 == mutual(i)
            plot([mapDataROI(idx2to1(idx1to2(i)),1) roverDataROI(idx1to2(i),1)],[mapDataROI(idx2to1(idx1to2(i)),2) roverDataROI(idx1to2(i),2)],'k','linewidth',2)
        end
    end
    
    
    % Find any data points that are missing from the reference (data set 1) in
    % the rover (data set 2)
    missing = zeros(length(roverDataROI),1);
    for i = 1:length(idx1to2)
        pos = find(i == idx2to1);
        if isempty(pos)
            missing(i) = 1;
            plot(mapDataROI(i,1),mapDataROI(i,2),'ko','Markersize',25)
        end
    end
    fprintf("%d features in the rover data set could not be matched to a corresponding map feature\n",sum(missing));
    
    % Find duplicate data points in the rover (data set 2) when compared with
    % the reference (data set 1)
    duplicate = zeros(length(roverDataROI),1);
    for i = 1:length(idx1to2)
        pos = find(i == idx2to1);
        if length(pos) > 1
            % Determine which of the duplicates is farther from the reference
            [~,minIdx] = min(sum((roverDataROI(pos,:)' - mapDataROI(i,:)').^2));
            for j = 1:length(pos)
                if j ~= minIdx
                    duplicate(pos(j)) = 1;
                    plot(roverDataROI(pos(j),1),roverDataROI(pos(j),2),'ks','Markersize',25)
                end
            end
        end
    end
    fprintf("%d features in the rover data set appear to be duplicates\n",sum(duplicate));
    pause(0.2);
end