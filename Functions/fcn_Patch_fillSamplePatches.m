function samplePatches = fcn_Patch_fillSamplePatches(varargin)
% fcn_Patch_fillSamplePatches
% Produces a pre-defined array of sample patche structures or loads an
% existing dataset from a file.
%
% FORMAT:
%
%       samplePatches = fcn_Patch_fillSamplePatches({filename})
%
% INPUTS:
%
%      (OPTIONAL INPUTS)
%      filename: a file name from which to load data
%
% OUTPUTS:
%
%      samplePatches: an cell array of datasets
%
% DEPENDENCIES:
%
%      (none)
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_Patch_fillSamplePatches.m for a full
%       test suite.
%
% This function was written on 2022_01_28 by C. Beal
% Questions or comments? cbeal@bucknell.edu

% Revision history:
%      2020_01_28
%      -- wrote the code

flag_do_debug = 0; % Flag to plot the results for debugging
flag_check_inputs = 1; % Flag to perform input checking
flag_load_from_file = 0; % Flag to determine whether to get data from a file

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
    if  nargin > 1
        error('Incorrect number of input arguments')
    end
end

if nargin > 0
    flag_load_from_file = 1;
    pathToDataFile = varargin{1};
    if ~exist(pathToDataFile,'file')
        error('No file exists at the specified location')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If the user has provided a file name
if flag_load_from_file
    s = load(pathToDataFile);
    sfields = fieldnames(s);
    if length(sfields) > 1
        error('Specified file is not formatted for loading a patch structure array');
    end
    samplePatches = getfield(s(1),sfields{1});
else
    % Create the empty scalar structure for the patches
    samplePatches = struct('id',{},'color',{},'primitive',{},'primparams',{},'aabb',{},'pointsX',{},'pointsY',{});
    
    % Create a test patch that is orange and very vaguely round
    samplePatches(1).id = 'bafded';       % Give it a unique identifier
    samplePatches(1).color = [1 0.6 0.3]; % Set the color to orange
    samplePatches(1).primitive = 'irregular'; % By default, set the primitive shape to 'irregular'
    % Add some points
    samplePatches(1).pointsX = [50.8640552995392; 58.2373271889401; 72.0622119815668; 68.6059907834101; 55.2419354838710];
    samplePatches(1).pointsY = [74.2335766423358; 66.9343065693431; 72.4817518248175; 85.3284671532847; 88.8321167883212];
    
    % Create a test patch that is gray and very vaguely rectangular
    samplePatches(2).id = 'cefdab';
    samplePatches(2).color = [0.6 0.6 0.6];
    samplePatches(2).primitive = 'irregular';
    % Add some points
    samplePatches(2).pointsX = [69.7580645161290; 86.3479262672811; 83.1221198156682; 79.6658986175115; 72.5230414746544; 69.2972350230415];
    samplePatches(2).pointsY = [8.83211678832116; 20.2189781021898; 28.1021897810219; 24.8905109489051; 21.3868613138686; 17.8832116788321];
end