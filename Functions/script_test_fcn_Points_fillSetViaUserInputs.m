% script_test_fcn_Points_fillSetViaUserInputs
% This is a script to exercise the function:
% fcn_Points_fillSetViaUserInputs.m
% This function was adapted on 2022_01_12 by C. Beal from S. Brennan's
% script_test_fcn_Path_fillPathViaUserInputs
% Questions or comments? sbrennan@psu.edu 

% Revisions
%     2022_01_12
%     -- adapted for XY datasets without path ordering

close all;
clc


%% BASIC example 1
if 1==0  % Intentionally comment this out so that doesn't autorun. Forces user to read the script!
    
    %% Code starts here
    fig_num = 1;
    h = figure(fig_num);
    hold on;
    
    num_iterations = input('How many paths do you want to draw? [Hit enter for default of 3]:','s');
    if isempty(num_iterations)
        num_iterations = 3;
    else
        num_iterations = str2double(num_iterations);
    end
    fprintf(1,'\n Filling in %.0d datasets.\n',num_iterations);
    fprintf(1,'Instructions: \n');
    fprintf(1,'Left click on the plot to create points. \n');
    fprintf(1,'Right click on the plot to remove points \n');
    fprintf(1,'Double click on the plot to end the dataset creation. \n');
    fprintf(1,'When the last dataset is completed, another plot will be created to show results. \n');
    
    
    % Initialize the paths_array
    clear dataset_array
    dataset_array{num_iterations} = [0 0];
    for i_set = 1:num_iterations
        
        % Set the title header
        UserData.title_header = sprintf('Dataset %.0d of %.0d',i_set,num_iterations);
        
        % Save the results
        set(gcf,'UserData',UserData);
        
        dataXY = fcn_Points_fillSetViaUserInputs(fig_num);
        dataset_array{i_set} = dataXY;
    end
    
    %% Plot the results

    fig_num = 13;
    title('Datasets recorded from user input')
    for i_set = 1:num_iterations
        plot(dataset_array{i_set},'.','Markersize',20);
    end
    
end