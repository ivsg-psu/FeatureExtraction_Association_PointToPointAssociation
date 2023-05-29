function datasets_out = fcn_Points_addRadialNoise(datasets,R,varargin)
% fcn_Points_addRadialNoise(datasets,R)
% Adds radial noise to the datasets using R (radius of max noise). The
% statistical distribution is sampled using a uniform sampling from 0 to R,
% and uniformly sampled along 0 to 2*pi.
%
% FORMAT: 
%
%       h = fcn_Points_addRadialNoise(datasets,R)
%
% INPUTS:
%
%      datasets: a cell array containing matricies of Nx2, of X and Y 
%      Note that each set will be modified separately.
%      R: radius of max noise
%
%      OPTIONAL INPUTS:
%      flag_noise_type: a flag to set noise types
%        1: uniform noise on R
%        2: random normal noise
%      fig_num: figure number
%
% OUTPUTS:
%
%      datasets_out: the input datasets will be modified and returned
%
% DEPENDENCIES:
%
%      fcn_Points_checkInputsToFunctions
%
% EXAMPLES:
%      
%       See the script: script_test_fcn_Points_addRadialNoise.m
%       for a full test suite. 
%
% This function was written on 2022_07_18 by Shashank Garikipati
% Questions or comments?

% Revision history:
% 2022_07_18
% -- wrote the code
% 2023_05_29 - sbrennan@psu.edu
% -- better comments, corrected the header
% -- vectorized code for speed
% -- added flag_noise_type

flag_do_debug = 0; % Flag to plot the results for debugging
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
    narginchk(2,4);
end

% Does user want to specify the distribution?
flag_noise_type = 1; % Flag to do a uniform distribution
if 3<= nargin
    temp = varargin{1};
    if ~isempty(temp)
        flag_noise_type = temp;
    end
end

% Does user want to show the plots?
flag_do_plot = 0; % Flag to do a plot
if 4 == nargin
    temp = varargin{end};
    if ~isempty(temp)
        fig_num = temp;
        figure(fig_num);
        flag_do_plot = 1;
    end
elseif flag_do_debug    
    fig = figure;
    fig_num = fig.Number;
    flag_do_plot = 1;
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
 
% Iterate through each of the datasets and add the radial noise


datasets_out = {zeros([size(datasets,2) 2])};
for i_set= 1:size(datasets,2)
    % How long is the current data set?
    N_data = length(datasets{i_set}(:,1));

   
    % Fill in X and Y values
    if flag_noise_type==1
        % Uniform random distribution on R
        rand_radius = R*rand(N_data,1);
        rand_angle = 2*pi*rand(N_data,1);
    elseif flag_noise_type==2
        % Random normal distribution on R
        rand_radius = R*randn(N_data,1);
        rand_angle = 2*pi*rand(N_data,1);
    else
        error('Unknown noise type specified.')
    end
    x = rand_radius.*cos(rand_angle) + datasets{i_set}(:,1);
    y = rand_radius.*sin(rand_angle) + datasets{i_set}(:,2);


    datasets_out{i_set} = [x y];
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
if flag_do_plot
    % Nothing in here yet
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file); 
end
end

