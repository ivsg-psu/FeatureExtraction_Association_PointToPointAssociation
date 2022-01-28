% script_test_fcn_Points_fillPointSetViaUserInputs
% This is a script to exercise the function:
% fcn_Path_fillPointSetViaUserInputs.m
% This function was adapted on 2022_01_28 by C. Beal from S. Brennan's
% script_test_fcn_Path_fillPathViaUserInputs.m
% Questions or comments? sbrennan@psu.edu 

% Revisions
%     2022_01_28
%     -- First write of the code

close all;
clc


%% BASIC example 1
if 1==1  % Intentionally comment this out so that doesn't autorun. Forces user to read the script!
    
    %% Code starts here
    fig_num = 1;
    h = figure(fig_num);
    hold on;
    
    num_iterations = input('How many point sets do you want to create? [Hit enter for default of 3]:','s');
    if isempty(num_iterations)
        num_iterations = 3;
    else
        num_iterations = str2double(num_iterations);
    end
    fprintf(1,'\n Filling in %.0d point sets.\n',num_iterations);
    fprintf(1,'Instructions: \n');
    fprintf(1,'Left click on the plot to create points. \n');
    fprintf(1,'Right click on the plot to remove points \n');
    fprintf(1,'Double click on the plot to end the set creation. \n');
    fprintf(1,'When the last set is completed, another plot will be created to show results. \n');
    
    
    % Initialize the paths_array
    clear point_array
    point_array{num_iterations} = [0 0];
    for i_set = 1:num_iterations
        
        % Set the title header
        UserData.title_header = sprintf('Set %.0d of %.0d',i_set,num_iterations);
        
        % Save the results
        set(gcf,'UserData',UserData);
        
        pointsXY = fcn_Points_fillPointSetViaUserInputs(fig_num);
        point_array{i_set} = pointsXY;
    end
    
    %% Plot the results
    clear data;
    
    % Plot the results
    fig_num = 13;
    fcn_Points_plotSetsXY(point_array);
    %fcn_Points_plotSetsXY(point_array,fig_num);
    
end