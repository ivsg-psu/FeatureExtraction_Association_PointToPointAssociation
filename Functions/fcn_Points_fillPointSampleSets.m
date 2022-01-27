function datasets_array = fcn_Points_fillSampleSets
% fcn_Points_fillSampleSets
% Produces dummy sample datasets. Note: can go into the function and change
% flag to allow user-selected datasets.
%
% FORMAT:
%
%       datasets_array = fcn_Points_fillSampleSets
%
% INPUTS:
%
%      (none)
%
% OUTPUTS:
%
%      datasets_array: an cell array of datasets
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_Points_fillSampleSets.m for a full
%       test suite.
%
% This function was adapted on 2022_01_12 by C. Beal from S. Brennan's
% script_test_fcn_Paths_fillSamplePaths
% Questions or comments? sbrennan@psu.edu

% Revision history:
%      2020_01_12 
%      -- adapted for XY datasets without path ordering

flag_do_debug = 0; % Flag to plot the results for debugging
flag_grab_user_inputs = 0; % Flag to allow user to click to draw paths
flag_check_inputs = 1; % Flag to perform input checking

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
end


%% check input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_check_inputs == 1
    % Are there the right number of inputs?
    if  nargin > 0
        error('Incorrect number of input arguments')
    end
    
end

%% Solve for the circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Build or load Set 1
if 1==flag_grab_user_inputs
    % Grab some points manually to create a starter dataset
    figure(1);
    axis([0 100 0 100]);
    
    % Get a set of starting points
    [X1,Y1] = ginput;
    datasets_array{1} = [X1 Y1];
else
    datasets_array{1} = [
       7.31566820276497,72.4817518248175;
       15.3801843317972,56.1313868613139;
       29.4354838709677,39.1970802919708;
       52.2465437788018,20.8029197080292;
       71.3709677419355,10.5839416058394;
       82.2004608294931,8.54014598540145];
end


% Build or load Set 2
if 1==flag_grab_user_inputs
    % Grab other points manually to create a test dataset
    figure(1);

    % Show prior results
    clf; hold on;
    for i_set = 1:length(datasets_array)
        plot(datasets_array{i_set}(:,1),datasets_array{i_set}(:,2),'-');
        text(datasets_array{i_set}(1,1),datasets_array{i_set}(1,2),'Start');
    end
    
    axis([0 100 0 100]);

    % Get a set of starting points
    [X2,Y2] = ginput;
    plot(X2,Y2)
    datasets_array{2} = [X2 Y2];
else
    datasets_array{2} = [
        13.7672811059908,75.9854014598540;
        20.9101382488479,62.2627737226277;
        34.7350230414747,44.4525547445255;
        52.7073732718894,30.1459854014599;
        75.7488479262673,12.6277372262774;
        84.9654377880184,8.54014598540145;
        90.7258064516129,7.08029197080292];
end

% % Build or load Set 3
% if 1==flag_grab_user_inputs
%     % Grab other points manually to create a test dataset
%     figure(1);
% 
%     % Show prior results
%     clf; hold on;
%     for i_set = 1:length(datasets_array)
%         plot(datasets_array{i_set}(:,1),datasets_array{i_set}(:,2),'-');
%         text(datasets_array{i_set}(1,1),datasets_array{i_set}(1,2),'Start');
%     end
%    
%     axis([0 100 0 100]);
%     
% 
%     % Get a set of starting points
%     [X3,Y3] = ginput;
%     plot(X3,Y3)
%     datasets_array{3} = [X3 Y3];
% else
%     datasets_array{3} = [
%         9.7926   10.6414
%         9.7926   15.8892
%         9.5622   19.0962
%         9.3318   19.9708
%         8.8710   21.7201
%         8.8710   22.8863
%         9.3318   24.0525
%         10.2535   25.2187
%         10.9447   26.6764
%         10.9447   26.9679
%         11.4055   28.4257
%         11.4055   29.8834
%         11.4055   31.0496
%         11.1751   31.6327
%         10.9447   33.6735
%         10.2535   35.7143
%         10.2535   36.8805
%         10.2535   38.3382
%         10.0230   40.3790
%         10.0230   41.5452
%         10.0230   43.2945
%         10.2535   45.9184
%         10.4839   49.4169
%         11.6359   52.9155
%         12.0968   55.5394
%         12.3272   58.4548
%         12.3272   60.4956
%         12.5576   64.2857
%         13.4793   69.2420
%         13.7097   71.8659
%         13.7097   74.4898
%         14.1705   76.8222
%         14.1705   77.9883
%         14.4009   79.7376
%         14.6313   82.3615
%         15.3226   83.8192
%         16.2442   85.2770
%         18.7788   88.1924
%         19.9309   89.3586
%         21.7742   89.9417
%         26.6129   89.9417
%         29.6083   90.2332
%         31.9124   90.5248
%         36.7512   91.1079
%         37.9032   91.1079
%         42.5115   91.1079
%         45.7373   91.1079
%         49.1935   91.6910
%         52.6498   91.6910
%         56.7972   90.5248
%         62.7880   90.5248
%         64.1705   89.6501
%         66.9355   85.5685
%         68.0876   82.6531
%         67.6267   80.3207
%         66.0138   77.1137
%         63.7097   75.0729
%         62.5576   73.6152
%         60.0230   69.5335
%         57.4885   65.4519
%         55.8756   62.8280
%         55.4147   59.6210
%         55.6452   58.4548
%         56.1060   57.5802
%         57.9493   56.1224
%         59.5622   55.5394
%         62.3272   56.4140
%         65.3226   57.8717
%         71.3134   57.8717
%         75.6912   57.5802
%         79.1475   56.4140
%         80.7604   55.8309
%         82.3733   54.9563
%         85.8295   45.6268
%         87.4424   44.1691
%         87.6728   40.6706
%         86.5207   39.5044
%         85.5991   32.2157
%         83.9862   28.7172
%         83.0645   26.6764
%         80.9908   24.0525
%         76.8433   21.7201
%         73.6175   20.5539
%         71.0829   19.0962
%         66.7051   18.8047
%         63.4793   19.6793
%         62.3272   20.8455
%         59.7926   22.8863
%         58.8710   24.6356
%         55.6452   26.9679
%         54.9539   28.7172
%         52.1889   31.9242
%         50.8065   33.3819
%         50.1152   33.6735
%         47.5806   34.5481
%         43.6636   34.5481
%         37.2120   35.1312
%         34.4470   35.1312
%         33.7558   37.4636
%         32.6037   41.2536
%         31.6820   43.8776
%         24.0783   47.3761
%         20.8525   47.9592
%         13.2488   48.5423
%         12.5576   49.1254
%         11.8664   52.9155
%         12.3272   57.2886
%         12.3272   59.9125
%         12.5576   62.2449
%         12.5576   65.4519
%         12.0968   69.2420
%         10.7143   75.6560
%         9.3318   78.5714
%         8.8710   82.6531
%         8.6406   84.4023
%         7.9493   86.4431
%         6.5668   88.4840
%         5.8756   89.6501];
% end


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
    % Prep a figure location
    close all
    figure(1);
    clf;
    hold on;
    grid minor;
    
    % Show result
    for i_set = 1:length(datasets_array)
        plot(datasets_array{i_set}(:,1),datasets_array{i_set}(:,2),'.','Markersize',20);
        text(datasets_array{i_set}(1,1),datasets_array{i_set}(1,2),'Start');
    end

end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end
end
