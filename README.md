# FeatureExtraction_Association_PointToPointAssociation

<!-- PROJECT LOGO -->
<br>
<p align="center">
  <h2 align="center">
    FeatureExtraction_Association_PointToPointAssociation
  </h2>
  <pre align="center">
        <img src=".\Images\PointToPointAssociation.jpg" alt="main mapgen picture" width="767" height="427">
        <!-- figcaption>Fig.1 - The typical progression of map generation.</figcaption -->
        <font size="-2">Photo by <a href="https://unsplash.com/@luism_arias?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">luis arias</a> on <a href="https://unsplash.com/photos/xyrz9dAGe6A?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a>
    </font>
  </pre>
</p>


<p align="left">
  MATLAB code implementation of a series of functions that associate spatial
  data. Specifically, functions are provided to determine matches between data
  sets of (X,Y) points, and store and compare groups of associated points.
</p>
<p align="center
">
    <br />
    <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation"><strong>Explore the docs »</strong></a>
    <br/>
    <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/tree/main/Documents">View Demo</a>
    ·
    <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/issues">Report Bug</a>
    ·
    <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/issues">Request Feature</a>
    <br />
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li>
      <a href="#structure">Repo Structure</a>
      <ul>
        <li><a href="#directories">Top-Level Directories</a></li>
        <li><a href="#dependencies">Dependencies</a></li>
      </ul>
    </li>
    <li>
    <a href="#functions">Functions</a>
    <ul>
        <li>
        <a href="#point-set-association-functions">Point-Set Association Functions</a>
         <ul>
         <li><a href="#fcn_points_plotsetsxy">fcn_Points_plotSetsXY</a></li>
          <li><a href="#fcn_points_fillpointsamplesets">fcn_Points_fillPointSampleSets</a></li>
          <li><a href="#fcn_points_fillpointsetviauserinputs">fcn_Points_fillPointSetViaUserInputs</a></li>
          <li><a href="#fcn_points_adjustpointsetstatistics">fcn_Points_adjustPointSetStatistics</a></li>
          <li><a href="#fcn_points_addradialnoise">fcn_Points_addRadialNoise</a></li>
          <li><a href="#fcn_points_calcpairstatistics">fcn_Points_calcPairStatistics</a></li>
          <li><a href="#fcn_points_pairxydata">fcn_Points_pairXYdata</a></li>
          <li><a href="#script_associatepointsincremental">script_AssociatePointsIncremental</a></li>
        </ul>
        </li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>


***
<!-- ABOUT THE PROJECT -->
## About The Project

<!--[![Product Name Screen Shot][product-screenshot]](https://example.com)-->

  MATLAB code implementation of a series of functions that associate spatial
  data. Specifically, functions are provided to determine matches between data
  sets of (X,Y) points, and store and compare groups of associated points.

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>

***

<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.

### Installation

1. Make sure to run MATLAB 2020b or higher. Why? The "digitspattern" command used in the DebugTools utilities was released late 2020 and this is used heavily in the Debug routines. If debugging is shut off, then earlier MATLAB versions will likely work, and this has been tested back to 2018 releases.

2. Clone the repo
   ```sh
   git clone https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation.git
   ```

3. Run the main code in the root of the folder (script_demo_PointToPointAssociation.m). This will download the required utilities for this code, unzip the zip files into a Utilities folder (.\Utilities), and update the MATLAB path to include the Utility locations. This install process will only occur the first time. Note: to force the install to occur again, delete the Utilities directory and clear all global variables in MATLAB (type: "clear global *").

4. Confirm it works! Run script_demo_PointToPointAssociation. If the code works, the script should run without errors. This script produces numerous example images such as those in this README file.

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>

***

<!-- STRUCTURE OF THE REPO -->
## Structure
### Directories
The following are the top level directories within the repository:
<ul>
	<li>Documents: Descriptions of the functionality and usage of the various MATLAB functions and scripts in the repository.</li>
	<li>Functions: The majority of the code for the point and patch association functionalities are implemented in this directory. All functions as well as test scripts are provided.</li>
	<li>Utilities: Dependencies that are utilized but not implemented in this repository are placed in the Utilities directory. These can be single files but are most often other cloned repositories.</li>
</ul>

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>

***

### Dependencies

* [Errata_Tutorials_DebugTools](https://github.com/ivsg-psu/Errata_Tutorials_DebugTools) - The DebugTools repo is used for the initial automated folder setup, and for input checking and general debugging calls within subfunctions. The repo can be found at: <https://github.com/ivsg-psu/Errata_Tutorials_DebugTools>

* [PathClass_v2023_02_01](https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary/blob/main/Releases/PathClass_v2023_02_01.zip?raw=true) - The PathClass repo is used for plotting lane markers and trajectory. The repo can be found at: <https://github.com/ivsg-psu/PathPlanning_PathTools_PathClassLibrary/blob/main/Releases/PathClass_v2023_02_01.zip?raw=true>

The dependencies are automatically installed by running the root master script

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>

***

<!-- FUNCTION DEFINITIONS -->
## Functions
The majority of the code for the point association functionalities are implemented in this directory. All functions as well as test scripts are provided.

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

### Point-Set Association Functions
<ul>
	<li>fcn_Points_fillPointSampleSets: a function to load some sample data sets to use for testing the other functions</li>
	<li>fcn_Points_fillPointSetViaUserInputs: a function that allows a user to create (X,Y) point sets by clicking in a figure with the mouse</li>
	<li>fcn_Points_plotSetsXY: a function that plots (X,Y) point sets with various options</li>
	<li>fcn_Points_pairXYdata: a function that associates the mutually closest points in two different point sets and returns the pairs as well as points which don't have an obvious mutual pair, in both directions</li>
	<li>fcn_Points_calcPairStatistics: a function that calculates the statistics for paired sets of points, returning RMS deviation, variance in point locations, and the offset between the centroids of the two point sets (a measurement of the systematic "shift" between two point sets)</li>
	<li>fcn_Points_adjustPointSetStatistics: a function to add 2D Gaussian noise and/or bias to a point set (e.g. to simulate sensor noise or bias) </li>
    <li>fcn_Points_addRadialNoise: Adds radial noise to the datasets using R (radius of max noise) </li>
</ul>

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>

***

#### fcn_Points_plotSetsXY

Plots the XY positions of all datasets existing in a data structure. This function is used to plot the data in the cell arrays.

FORMAT: 
```MATLAB
      h = fcn_Points_plotSetsXY(datasets,{fig_num})
```
**INPUTS:**

datasets: a structure array containing subfields of X and Y 
coordinates in the following form:

    datasets{i_set}.X   
    datasets{i_set}.Y

fig_num is an optional input.

Note that i_set denotes a data set structure. Each set can be
plotted separately and all the dataset can also be plotted in the same figure.  

**OUTPUTS:**

h: a handle to the resulting figure

**Examples:**

```MATLAB
% Example 1
% Plots XY positions of one dataset
fig_num = 8;

dataset = {};

dataset{1} = [1,21;
              2,2;
              3,3;
              4,1;
              5,4;
              6,32;
              7,0.5];

fcn_Points_plotSetsXY(dataset,fig_num)
title('fcn_Points_plotSetsXY: Example plotting XY positions of one dataset','Interpreter','none')
```
<pre align="center">
  <img src=".\Images\fcn_Points_plotSetsXY_Ex1.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Example plotting XY positions of one dataset.</figcaption>
</pre>

```MATLAB
% Example 2
% Plots XY positions of three datasets

fig_num = 9;

dataset2 = {};
dataset2{1} = [1,21;
              2,2;
              3,3;
              4,1;
              5,4;
              6,32;
              7,0.5];

dataset2{2} = dataset{1}*1.8+0.3;

dataset2{3} = [9,2;
              8,0;
              2,5;
              1,1;
              6,4;
              10,31;
              2,1];

fcn_Points_plotSetsXY(dataset2,fig_num)
title('fcn_Points_plotSetsXY: Example plotting XY positions of three datasets','Interpreter','none')
```
<pre align="center">
  <img src=".\Images\fcn_Points_plotSetsXY_Ex2.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Example plotting XY positions of three datasets.</figcaption>
</pre>

```MATLAB
% Example 3
% Plots from a sample data set inside fcn_Points_fillPointSampleSets

fig_num = 10;
datasets3 = fcn_Points_fillPointSampleSets;

fcn_Points_plotSetsXY(datasets3,fig_num);
title('fcn_Points_plotSetsXY: Example plotting sample data set inside fcn_Points_fillPointSampleSets','Interpreter','none')
for i_set = 1:length(datasets3)
    text(datasets3{i_set}(1,1),datasets3{i_set}(1,2),'Start');
end

```
<pre align="center">
  <img src=".\Images\fcn_Points_plotSetsXY_Ex3.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption> Example plotting sample data set in fcn_Points_fillPointSampleSets.</figcaption>
</pre>

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### fcn_Points_fillPointSampleSets

Produces pre-defined sample datasets or loads an existing dataset from a
file.

FORMAT:

```MATLAB
      datasets_array = fcn_Points_fillPointSampleSets({filename})
```
INPUTS:

(OPTIONAL INPUTS)
filename: a file name from which to load data

OUTPUTS:

datasets_array: an cell array of datasets

**Examples:** 

```MATLAB
% Setup
% Call the function to fill in an array of "dataset" type
datasets_array = fcn_Points_fillPointSampleSets;


% Example 1
% We can save one of these as a single "dataset"
fig_num = 1;
single_path = {datasets_array{1}};

fcn_Points_plotSetsXY(single_path,fig_num);
text(single_path{1}(1,1),single_path{1}(1,2),'Start');
title('fcn_Points_plotSetsXY: Example plotting of one dataset ','Interpreter','none')  
```
<pre align="center">
  <img src=".\Images\fcn_Points_fillPointSampleSets_Ex1.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Plot for the 1st dataset in fcn_Points_fillPointSampleSets.</figcaption>
</pre>

```MATLAB

% Example 2
% illustrating the different XY points on different plots

for i_set = 1:length(datasets_array)
    fig_num = i_set+1;
    fcn_Points_plotSetsXY(datasets_array(i_set),fig_num);
    text(datasets_array{i_set}(1,1),datasets_array{i_set}(1,2),'Start');
    title(sprintf('fcn_Points_plotSetsXY: Example plotting of dataset: %.0d',i_set),'Interpreter','none'); 
end

```
<pre align="center">
  <img src=".\Images\fcn_Points_fillPointSampleSets_Ex2_1.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Plot for the 1st dataset in fcn_Points_fillPointSampleSets.</figcaption>
  
</pre>

<pre align="center">
  <img src=".\Images\fcn_Points_fillPointSampleSets_Ex2_2.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Plot for the 2nd dataset in fcn_Points_fillPointSampleSets.</figcaption>
</pre>

<pre align="center">
  <img src=".\Images\fcn_Points_fillPointSampleSets_Ex2_3.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Plot for the 3rd dataset in fcn_Points_fillPointSampleSets.</figcaption>
</pre>

```MATLAB
fig_num = 5;
fcn_Points_plotSetsXY(datasets_array,fig_num);
title('fcn_Points_plotSetsXY: Example plotting of all datasets in one call','Interpreter','none')


```

<pre align="center">
  <img src=".\Images\fcn_Points_fillPointSampleSets_Ex3.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Plot for the all the datasets in fcn_Points_fillPointSampleSets.</figcaption>
</pre>

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### fcn_Points_fillPointSetViaUserInputs

A function for the user to click on the figure to generate XY data.
Points are collected and plotted until the user double clicks. If the
user right-clicks anywhere in the plot, the last point is deleted. Once
the user double-clicks, the results are output from the function.

FORMAT:

```MATLAB
     dataXY = fcn_Points_fillPointSetViaUserInputs({fig_num})
```

**INPUTS:**

(OPTIONAL INPUTS)
fig_num: an integer specifying which figure to use

**OUTPUTS:**

dataXY: matrix (Nx2) representing the X and Y points that the user
clicked on the map

**Examples:** 

```MATLAB
% Shows step by step instructions about using function

fig_num = 6;
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
```
<pre align="center">
  <img src=".\Images\fcn_Points_fillPointSetViaUserInputs_Ex1_1.png" alt="fcn_Points_fillPointSampleSets picture" width="500" height="100">
  <figcaption>Instructions and selecting number of point sets to create.</figcaption>
</pre>

```MATLAB
    
% Initialize the paths_array
clear point_array
point_array{num_iterations} = [0 0];
for i_set = 1:num_iterations
    
    % Set the title header
    UserData.title_header = sprintf('Set %.0d of %.0d',i_set,num_iterations);
    
    % Save the results
    set(gcf,'UserData',UserData);
    
    pointsXY = fcn_Points_fillPointSetViaUserInputs(fig_num);
    title('fcn_Points_fillPointSetViaUserInputs: Plotting the points via user inputs','Interpreter','none')
    point_array{i_set} = pointsXY;
end
```
<pre align="center">
  <img src=".\Images\fcn_Points_fillPointSetViaUserInputs_Ex1_2.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Plotting the set1 points via user inputs.</figcaption>
</pre>

<pre align="center">
  <img src=".\Images\fcn_Points_fillPointSetViaUserInputs_Ex1_3.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Plotting the set2 points via user inputs.</figcaption>
</pre>

```MATLAB
clear data, close;
% Plot the results
fig_num = 7;
fcn_Points_plotSetsXY(point_array,fig_num);
title('Plotting the point array (points plotted via user inputs) in a new figure','Interpreter','none')

```
<pre align="center">
  <img src=".\Images\fcn_Points_fillPointSetViaUserInputs_Ex2.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Plotting the point array (points plotted via user inputs).</figcaption>
</pre>

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### fcn_Points_adjustPointSetStatistics

Adjusts the XY positions of all datasets existing in a data structure by
the provided values of bias, noise mean and noise variance.

FORMAT: 

```MATLAB
      h = fcn_Points_adjustPointSetStatistics(datasets,bias,noiseMean,noiseVariance)
```
**INPUTS:**

datasets: a structure array containing subfields of X and Y 
coordinates in the following form:
<ul>
  datasets{i_set}.X
    
  datasets{i_set}.Y
</ul>
Note that i_set addresses a particular data set structure. Each set
will be modified separately.

<ul>

bias: a matrix that defines constant X and Y offsets for each data 
        set, hence an N x 2 matrix where N is the number of data sets

noiseMean: a matrix that defines the X and Y means for the noise
          to be added to each of the each data sets, hence N x 2

noiseVariance: a matrix that defines the X and Y variance for the
              noise to be added to each of the each data sets, hence N x 2
</ul>

**OUTPUTS:**

datasets: the input datasets will be modified and returned

**Examples:** 

```MATLAB
% Load some test data sets
origXYdatasets = fcn_Points_fillPointSampleSets;
biasx = -2.2;
biasy = 1.2;

noiseMean = zeros(1,2);
noiseVariance = zeros(1,2);
%Corrupt one set of test data with only a systematic bias (for each)
biasedXYdataset = fcn_Points_adjustPointSetStatistics(origXYdatasets(1),[biasx biasy],noiseMean,noiseVariance);

% Plotting the original data and biased data
% fh(1) = figure;
fig_num = 12;
figure(fig_num)
clf;
axis equal
grid on
xlabel('x [m]')
ylabel('y [m]')
    
horig(1) = fcn_Points_plotSetsXY(origXYdatasets(1),figure(fig_num))
set(horig(1),'color','blue');
hcorr(1) = fcn_Points_plotSetsXY(biasedXYdataset,figure(fig_num))
set(hcorr(1),'marker','*')
set(hcorr(1),'color','blue');
legend([horig(1) hcorr(1)],{'Original data','Biased data'})
title('fcn_Points_adjustPointSetStatistics: Plot of the original and biased datasets','Interpreter','none');
```
<pre align="center">
  <img src=".\Images\fcn_Points_adjustPointSetStatistics.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Plot of the original and biased datasets.</figcaption>
</pre>

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### fcn_Points_addRadialNoise

Adds radial noise to the datasets using R (radius of max noise)

FORMAT: 
```MATLAB
      h = fcn_Points_addRadialNoise(datasets,bias,R)
```
**INPUTS:**

datasets: a structure array containing subfields of X and Y 
coordinates in the following form:
<ul>
  datasets{i_set}.X

  datasets{i_set}.Y
</ul>
Note that i_set addresses a particular data set structure. Each set
will be modified separately.
R: radius of max noise


OPTIONAL INPUTS:
fig_num: figure number

**OUTPUTS:**

datasets_out: the input datasets will be modified and returned

**Examples:**
```MATLAB
% Load some test data sets
origXYdatasets = fcn_Points_fillPointSampleSets;

% Example 1
% shows addition of radial noise

fig_num = 14;
new_set1 = {[0.0 0.0]};


RadiusMaxNoise = 0.5;
for i =  1:100
    datasets_out = fcn_Points_addRadialNoise(new_set1,RadiusMaxNoise);
    fcn_Points_plotSetsXY(new_set1,fig_num);
    fcn_Points_plotSetsXY(datasets_out,fig_num);
    title('fcn_Points_addRadialNoise: Demonstrates the addition of radial noise','Interpreter','none');
    axis([-1 1 -1 1])
    axis square
end
```
<pre align="center">
  <img src=".\Images\fcn_Points_addRadialNoise_Ex1.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Demonstrates the addition of radial noise.</figcaption>
</pre>

```MATLAB
% Example 2
% Addition of radial noise for one data set

fig_num = 15;
RadiusMaxNoise = 1;
new_set2 = {origXYdatasets{1}}; %#ok<CCAT1> 
datasets_out2 = fcn_Points_addRadialNoise(new_set2,RadiusMaxNoise);
fcn_Points_plotSetsXY(new_set2,fig_num); 
fcn_Points_plotSetsXY(datasets_out2,fig_num);
title('fcn_Points_addRadialNoise: Adding noise to the one of the datasets','Interpreter','none')
```
<pre align="center">
  <img src=".\Images\fcn_Points_addRadialNoise_Ex2.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Adding noise to the one of the datasets.</figcaption>
</pre>

```MATLAB
% Example 3
% Addition of radial noise for one data set
fig_num = 16;
RadiusMaxNoise = 1;

new_set3 = origXYdatasets;
datasets_out3 = fcn_Points_addRadialNoise(new_set3,RadiusMaxNoise);
fcn_Points_plotSetsXY(new_set3,fig_num);
fcn_Points_plotSetsXY(datasets_out3,fig_num);
title('fcn_Points_addRadialNoise: Adding noise to all the datasets','Interpreter','none')

```
<pre align="center">
  <img src=".\Images\fcn_Points_addRadialNoise_Ex3.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Adding noise to all the datasets.</figcaption>
</pre>

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>

***

#### fcn_Points_calcPairStatistics

Calculates the statistics of the errors between two sets of paired XY
data points and an overall X,Y map shift to match the centroids of the
data sets.

FORMAT:

```MATLAB
      [errRMS,errVar,meanShift] = fcn_Points_calcPairStatistics(pairedXYdata)
```
**INPUTS:**

pairedXYdata: an N x 4 matrix with [X1 Y1 X2 Y2] data in each row.

**OUTPUTS:**
<ul>

  errRMS: the RMS distances associated with the point-to-point match errors.     
  

  errVar: the variance in the point-to-point match error distances.
 
  meanShift: the [X Y] shift that would cause the centroids of the
            two data sets to match.
</ul>

**Examples**

```MATLAB
% Plotting a sine wave
t=0:0.01:1; 
f=1; 
x=sin(2*pi*f*t);

figure(88); 
plot(t,x,'LineWidth',1.5);
title('Sinusoidal wave with a frequency of 1 Hz','Interpreter','none');
```
<pre align="center">
  <img src=".\Images\fcn_Points_calcPairStatistics_Ex1.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Sinusoidal wave with a frequency of 1 Hz.</figcaption>
</pre>

```MATLAB

% creating a dataset using t and x values 
data_set = [t;x]';

% chainging the data set into cell array
XYdata = {data_set};

% Adding bias to XYdata
biasx = -2.3;
biasy = 0.9;

noiseMean = zeros(1,2);
noiseVariance = zeros(1,2);

biasedXYdata = fcn_Points_adjustPointSetStatistics(XYdata,[biasx biasy],noiseMean,noiseVariance);


% Plotting the original and biased data on the same figure
fig_num = 21;
figure(fig_num)
axis equal
grid on
xlabel('t')
ylabel('x')

Orig = fcn_Points_plotSetsXY(XYdata,figure(fig_num));
set(Orig,'color','blue');
biaS = fcn_Points_plotSetsXY(biasedXYdata,figure(fig_num));
set(biaS,'marker','*')
set(biaS,'color','blue');

legend([Orig biaS],{'Original data','Biased data'})
title('fcn_Points_adjustPointSetStatistics: Plot of the original and biased datasets','Interpreter','none');

% Run the statistics on the original and biased data

OrignBiasData = [XYdata{1},biasedXYdata{1}];

[errRMS,errVar,meanShift] = fcn_Points_calcPairStatistics(OrignBiasData);

fprintf(1,'\n The RMS error is %.4f. \n',errRMS);
fprintf(1,'The variance is %.4f. \n',errVar);
fprintf(1,'The meanShift is %.4f. \n',meanShift);

```
<pre align="center">
  <img src=".\Images\fcn_Points_calcPairStatistics_Ex2.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Plot of the original and biased datasets.</figcaption>
</pre>

```MATLAB

% Noisy data
RadiusMaxNoise = 0.5;
NoisyXYData = fcn_Points_addRadialNoise(XYdata,RadiusMaxNoise);

axis equal
grid on
fig_num = 51;    
Orig = fcn_Points_plotSetsXY(XYdata,fig_num);
set(Orig,'color','blue');
noisY = fcn_Points_plotSetsXY(NoisyXYData,fig_num);
set(noisY,'marker','*')
set(noisY,'color','red')

legend([Orig noisY],{'Original data','Noisy data'})
title('fcn_Points_addRadialNoise: Plot of the original and noisy datasets','Interpreter','none');

OrignNoisyData = [XYdata{1},NoisyXYData{1}];

[errRMS,errVar,meanShift] = fcn_Points_calcPairStatistics(OrignNoisyData);

fprintf(1,'\n The RMS error is %.4f. \n',errRMS);
fprintf(1,'The variance is %.4f. \n',errVar);
fprintf(1,'The meanShift is %.4f. \n',meanShift);

```
<pre align="center">
  <img src=".\Images\fcn_Points_calcPairStatistics_Ex3.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Plot of the original and noisy datasets.</figcaption>
</pre>

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>

***

#### fcn_Points_pairXYdata

Determines the closest pairs of matching X,Y data points from two
different data sets of (potentially) different length. The algorithm uses
a "mutual" pairing approach, where the nearest neighbor point is computed
for each point in each data set and then the matches between the closest
points are paired up. The matched pairs are listed in a single matrix,
followed by unmatched pairs from data set A and then data set B. The
unmatched data points are padded with NaN values. The number of matched
points as well as number of unmatched points from data set A and data set
B are returned.

FORMAT:
```MATLAB
      [pairedXYdata, numMatches, nonMatchesA, nonMatchesB] = fcn_Points_pairXYdata(xyDataA,xyDataB,(maxDist))
```

**INPUTS:**

xyDataA: an N x 2 matrix with [X Y] data in each row.

xyDataB: an N x 2 matrix with [X Y] data in each row.

(OPTIONAL INPUTS)
maxDist: a maximum radius outside of which to reject matched pairs

**OUTPUTS:**

pairedXYdata: an N x 4 matrix with [X1 Y1 X2 Y2] data in each row

numMatches: the number of mutual matches in data sets A and B

nonMatchesA: the number of points in data set A without a mutual match

nonMatchesB: the number of points in data set B without a mutual match

**Examples**
```MATLAB
% Load up some data (simple xy points for now)
load testDatasetVehicle.mat

% Call the pairing function to obtain the matrix of paired XY data without
% a limiting radius
%[pairedXYdata, numMatches, nonMatchesA, nonMatchesB] = fcn_Points_pairXYdata(xyData{1},xyData{2});

% Define a maximum radius within which to consider data points as pairs
pairRadius = 1.0;

% Call the pairing function to obtain the matrix of paired XY data
[pairedXYdata, numMatches, nonMatchesA, nonMatchesB] = fcn_Points_pairXYdata(xyData{1},xyData{2},pairRadius);


% Now plot the matches with stars in various colors (to identify the pairs)
figure(11);
clf
hold on
quiver(pairedXYdata(:,1),pairedXYdata(:,2),pairedXYdata(:,3)-pairedXYdata(:,1),pairedXYdata(:,4)-pairedXYdata(:,2),0,'k','linewidth',1)

for i = 1:numMatches
    plot([pairedXYdata(i,1) pairedXYdata(i,3)],[pairedXYdata(i,2) pairedXYdata(i,4)],'*')
end
% Also plot the non-matches in data set A with circles in various colors
for i = 1+numMatches:nonMatchesA+numMatches
    plot(pairedXYdata(i,1),pairedXYdata(i,2),'o')
end
% Lastly, plot the non-matches in data set B with squares in various colors
for i = 1+numMatches+nonMatchesA:nonMatchesA+numMatches+nonMatchesA
    plot(pairedXYdata(i,3),pairedXYdata(i,4),'d')
end

title('fcn_Points_pairXYdata: Pairing the nearest points','Interpreter','none')
```
<pre align="center">
  <img src=".\Images\fcn_Points_pairXYdata_Ex1.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Pairing the nearest points.</figcaption>
</pre>

```MATLAB
% Load up some data (simple xy points for now)
load testDatasetVehicle.mat

% Call the pairing function to obtain the matrix of paired XY data without
% a limiting radius
%[pairedXYdata, numMatches, nonMatchesA, nonMatchesB] = fcn_Points_pairXYdata(xyData{1},xyData{2});

% Define a maximum radius within which to consider data points as pairs
pairRadius = 3.0;

% Call the pairing function to obtain the matrix of paired XY data
[pairedXYdata, numMatches, nonMatchesA, nonMatchesB] = fcn_Points_pairXYdata(xyData{1},xyData{2},pairRadius);


% Now plot the matches with stars in various colors (to identify the pairs)
figure(11);
clf
hold on
quiver(pairedXYdata(:,1),pairedXYdata(:,2),pairedXYdata(:,3)-pairedXYdata(:,1),pairedXYdata(:,4)-pairedXYdata(:,2),0,'k','linewidth',1)

for i = 1:numMatches
    plot([pairedXYdata(i,1) pairedXYdata(i,3)],[pairedXYdata(i,2) pairedXYdata(i,4)],'*')
end
% Also plot the non-matches in data set A with circles in various colors
for i = 1+numMatches:nonMatchesA+numMatches
    plot(pairedXYdata(i,1),pairedXYdata(i,2),'o')
end
% Lastly, plot the non-matches in data set B with squares in various colors
for i = 1+numMatches+nonMatchesA:nonMatchesA+numMatches+nonMatchesA
    plot(pairedXYdata(i,3),pairedXYdata(i,4),'d')
end

title('fcn_Points_pairXYdata: Pairing the nearest points for larger pairRadius','Interpreter','none')
```
<pre align="center">
  <img src=".\Images\fcn_Points_pairXYdata_Ex2.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Pairing the nearest points for larger pairRadius.</figcaption>
</pre>

```MATLAB
CalcpairedXYdata = pairedXYdata(1:numMatches,:);

% Run the statistics on the paired data
[errRMS,errVar,meanShift] = fcn_Points_calcPairStatistics(CalcpairedXYdata);

fprintf(1,'\n The RMS error is %.4f. \n',errRMS);
fprintf(1,'The variance is %.4f. \n',errVar);
fprintf(1,'The meanShift is %.4f. \n',meanShift);


% Post-process the shift to extract a distance and an angle for easy visual
% checking of the data
shiftDist = norm(meanShift,2);
shiftAngle = atan2(meanShift(2),meanShift(1));

% Load some test data sets
origXYdatasets = fcn_Points_fillPointSampleSets;
biasx = -1.2;
biasy = 0.6;

noiseMean = zeros(1,2);
noiseVariance = zeros(1,2);
%Corrupt one set of test data with only a systematic bias (for each)
biasedXYdataset = fcn_Points_adjustPointSetStatistics(origXYdatasets(1),[biasx biasy],noiseMean,noiseVariance);


% Call the pairing function to obtain pair the original data with the
% biased data
[pairedXYdataBias, numMatchesBias, nonMatchesABias, nonMatchesBBias] = fcn_Points_pairXYdata(origXYdatasets{1},biasedXYdataset{:});

% Calculate the statistics for the biased data set relative to the original
[errRMSBias,errVarBias,meanShiftBias] = fcn_Points_calcPairStatistics(pairedXYdataBias(1:numMatchesBias,:));

% fh(1) = figure;
fig_num = 12;
figure(fig_num)
set(gcf,'Position',[100 100 1200 500])
clf;
% Plot to provide a visual inspection of the bias
subplot(1,2,1)

grid on
xlabel('x [m]')
ylabel('y [m]')
    
horig(1) = fcn_Points_plotSetsXY(origXYdatasets(1),figure(fig_num));
set(horig(1),'color','blue');
hcorr(1) = fcn_Points_plotSetsXY(biasedXYdataset,figure(fig_num));
set(hcorr(1),'marker','*')
set(hcorr(1),'color','blue');
legend([horig(1) hcorr(1)],{'Original data','Biased data'})

subplot(1,2,2)
% Plotting original data paired with biased data
% figure(13);
% clf

grid on
hold on
xlabel('x [m]')
ylabel('y [m]')
quiver(pairedXYdataBias(:,1),pairedXYdataBias(:,2),pairedXYdataBias(:,3)-pairedXYdataBias(:,1),pairedXYdataBias(:,4)-pairedXYdataBias(:,2),0,'k','linewidth',1)
for i = 1:size(pairedXYdataBias,1)
    plot([pairedXYdataBias(i,1) pairedXYdataBias(i,3)],[pairedXYdataBias(i,2) pairedXYdataBias(i,4)],'*')
end

sgtitle('fcn_Points_pairXYdata with Bias: Plotting the biased dataset and pairing the biased datatset with original','Interpreter','none')
```
<pre align="center">
  <img src=".\Images\fcn_Points_pairXYdata_Ex3.jpg" alt="fcn_Points_fillPointSampleSets picture" width="600" height="300">
  <figcaption>Plotting the biased dataset and pairing the biased datatset with original.</figcaption>
</pre>


<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### script_AssociatePointsIncremental

This script simulates a vehicle following a trajectory by detecting the nearest points usings KNN algoorithm.

```MATLAB
clear all
close all

% define the dot product
dot = @(x,y) x(:)'*y(:);

% Assign a region of interest around the rover
ROIradius = 10;

% Assign a radial standard deviation value
sigma = 3; % meters

load testDatasetVehicle.mat

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
```
<pre align="center">
  <img src=".\Images\script_AssociatePointsIncremental_Ex1.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Plot for all the datasets in testDatasetVehicle.mat.</figcaption>
</pre>

```MATLAB
for locInd = 1:length(xyData{3})
    % Set up for another round of drawing on the figure
    figure(8)
    clf
    hold on
    axis ([0 100 0 100])
    
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
    pause(0.3);
 end
```

<pre align="center">
  <img src=".\Images\script_AssociatePointsIncremental_Ex1.jpg" alt="fcn_Points_fillPointSampleSets picture" width="600" height="300">
  <figcaption>Plotting the biased dataset and pairing the biased datatset with original.</figcaption>
</pre>

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***


<!-- USAGE EXAMPLES -->
## Usage
<!-- Use this space to show useful examples of how a project can be used.
Additional screenshots, code examples and demos work well in this space. You may
also link to more resources. -->

1. Open MATLAB and navigate to the Functions directory

2. Run any of the various test scripts, such as
   ```sh
   script_test_fcn_Points_pairXYdata
   ```
   or
   ```sh
   script_test_fcn_Patch_checkCollisions
   ```
_For more examples, please refer to the [Documentation] 

https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/tree/main/Documents)_

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

<!-- CONTACT -->
## Contact

Prof. Sean Brennan - sbrennan@psu.edu

Project Link: [FeatureExtraction_Association_PointToPointAssociation](https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation)

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation.svg?style=for-the-badge
[contributors-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation.svg?style=for-the-badge
[forks-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/network/members
[stars-shield]: https://img.shields.io/github/stars/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation.svg?style=for-the-badge
[stars-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/stargazers
[issues-shield]: https://img.shields.io/github/issues/ivsg-psu/reFeatureExtraction_Association_PointToPointAssociationpo.svg?style=for-the-badge
[issues-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/issues
[license-shield]: https://img.shields.io/github/license/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation.svg?style=for-the-badge
[license-url]: https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/blob/master/LICENSE.txt

