function datasets = fcn_Points_adjustPointSetStatistics(datasets,bias,noiseMean,noiseVariance)
% fcn_Points_adjustPointSetStatistics
% Adjusts the XY positions of all datasets existing in a data structure by
% the provided values
%
% FORMAT: 
%
%       h = fcn_Points_adjustPointSetStatistics(datasets,bias,noiseMean,noiseVariance)
%
% INPUTS:
%
%      datasets: a structure array containing subfields of X and Y 
%      coordinates in the following form:
%           datasets{i_set}.X
%           datasets{i_set}.Y
%      Note that i_set addresses a particular data set structure. Each set
%      will be modified separately.
%      bias: a matrix that defines constant X and Y offsets for each data 
%            set, hence an N x 2 matrix where N is the number of data sets
%      noiseMean: a matrix that defines the X and Y means for the noise
%            to be added to each of the each data sets, hence N x 2
%      noiseVariance: a matrix that defines the X and Y variance for the
%            noise to be added to each of the each data sets, hence N x 2
%
% OUTPUTS:
%
%      datasets: the input datasets will be modified and returned
%
% DEPENDENCIES:
%
%      fcn_Points_checkInputsToFunctions
%
% EXAMPLES:
%      
%       See the script: script_test_fcn_Points_adjustPointSetStatistics.m
%       for a full test suite. 
%
% This function was written on 2022_01_12 by C. Beal
% Questions or comments? cbeal@bucknell.edu 

% Revision history:
%     2022_01_28
%     -- wrote the code

flag_do_debug = 0; % Flag to plot the results for debugging
flag_this_is_a_new_figure = 1; % Flag to check to see if this is a new figure
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
NumSets = length(datasets);

if flag_check_inputs == 1
    % Are there the right number of inputs?
    if 4 ~= nargin
        error('Incorrect number of input arguments')
    end
    
    % Check the data input for size consistency
    if NumSets ~= size(bias,1)
        error('Bias vector should have the same number of rows as the number of data sets.')
    end
    if 2 ~= size(bias,2)
        error('Bias vector should have two columns for the X and Y bias for each data set.')
    end
    if NumSets ~= size(noiseMean,1)
        error('Noise mean vector should have the same number of rows as the number of data sets.')
    end
    if 2 ~= size(noiseMean,2)
        error('Noise mean vector should have two columns for the X and Y bias for each data set.')
    end
    if NumSets ~= size(noiseVariance,1)
        error('Noise variance vector should have the same number of rows as the number of data sets.')
    end
    if 2 ~= size(noiseVariance,2)
        error('Noise variance vector should have two columns for the X and Y bias for each data set.')
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
 
% Iterate through each of the datasets and add the bias and noise
for i_set= 1:NumSets
    % Generate the appropriate noise
    noiseMatrix = mvnrnd(noiseMean(i_set,:),diag(noiseVariance(i_set,:)),size(datasets{i_set},1));
    % Add the bias to the data
    datasets{i_set} = [datasets{i_set}(:,1) + bias(i_set,1), datasets{i_set}(:,2) + bias(i_set,2)];
    % Add the noise to the data
    datasets{i_set} = datasets{i_set} + noiseMatrix;
end


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

