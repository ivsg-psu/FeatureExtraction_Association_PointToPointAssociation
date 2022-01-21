% Script to process data sets and test for congruency

clear all
close all

% Generate two data sets via interactive plot
% figure(9);
% hold on;
% axis equal
% % Set the title header
% UserData.title_header = sprintf('Dataset %.0d of %.0d',1,2);
% % Save the results
% set(gcf,'UserData',UserData);
% inputData = fcn_Dataset_fillSetViaUserInputs(9);
% xyData{1} = inputData;
% % Set the title header
% UserData.title_header = sprintf('Dataset %.0d of %.0d',2,2);
% % Save the results
% set(gcf,'UserData',UserData);
% inputData = fcn_Dataset_fillSetViaUserInputs(9);
% xyData{2} = inputData;
% 
% % Assign a radial standard deviation value
% sigma = 3; % meters
% 
% % Load up some data (simple xy points for now)
% %xyData = fcn_Dataset_fillSampleSets;
load testDataset1.mat
% 
%% Plot the data (known to be two sets for now)
figure(1)
clf
hold on
for i = 1:length(xyData{1})
    plot(xyData{1}(i,1),xyData{1}(i,2),'b.','Markersize',20);
    text(xyData{1}(i,1)+2,xyData{1}(i,2),sprintf('%d',i))
end
for i = 1:length(xyData{2})
    plot(xyData{2}(i,1),xyData{2}(i,2),'r.','Markersize',20);
    text(xyData{2}(i,1)+2,xyData{2}(i,2),sprintf('%d',i))
end


% Binary compatibility analysis (mutually nearest neighbor points)
% First, associate data points in the sets using a nearest neighbor search
[idx1to2,dist1] = knnsearch(xyData{2},xyData{1});
[idx2to1,dist2] = knnsearch(xyData{1},xyData{2});

% Find mutual matches between nearest neighbors. These points are "most" in
% agreement in that they are mutual nearest neighbor matches
mutual = nan(length(idx1to2),1);
for i = 1:length(idx1to2)
    if(idx2to1(idx1to2(i)) == i)
        mutual(i) = 1;
    else
        mutual(i) = 0;
    end
end
% Indicate on the plot that these points match
for i = 1:length(mutual)
    if 1 == mutual(i)
        plot([xyData{1}(idx2to1(idx1to2(i)),1) xyData{2}(idx1to2(i),1)],[xyData{1}(idx2to1(idx1to2(i)),2) xyData{2}(idx1to2(i),2)],'k','linewidth',2)
    end
end

% Determine the binary compatibility of the data sets (percentage of points
% that are mutually nearest neighbors)
binMetric = sum(mutual)/length(mutual);
fprintf("%.2f percent of map features matched in rover data set\n",100*binMetric); 
% Determine the rms error between mutually compatible points
rmsError = sqrt(mean(dist1(mutual == 1).^2));
fprintf("RMS error for matched features is %.2f units\n",rmsError);

% Find any data points that are missing from the reference (data set 1) in
% the rover (data set 2)
missing = zeros(length(xyData{2}),1);
for i = 1:length(idx1to2)
    pos = find(i == idx2to1);
    if isempty(pos)
        missing(i) = 1;
        plot(xyData{1}(i,1),xyData{1}(i,2),'ko','Markersize',25)
    end
end
fprintf("%d features in the rover data set could not be matched to a corresponding map feature\n",sum(missing));

% Find duplicate data points in the rover (data set 2) when compared with
% the reference (data set 1)
duplicate = zeros(length(xyData{2}),1);
for i = 1:length(idx1to2)
    pos = find(i == idx2to1);
    if length(pos) > 1
        % Determine which of the duplicates is farther from the reference
        [~,minIdx] = min(sum((xyData{2}(pos,:)' - xyData{1}(i,:)').^2));
        for j = 1:length(pos)
            if j ~= minIdx
                duplicate(pos(j)) = 1;
                %plot(xyData{2}(pos(j),1),xyData{2}(pos(j),2),'ks','Markersize',25)
            end
        end
    end
end
fprintf("%d features in the rover data set appear to be duplicates\n",sum(duplicate));