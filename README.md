# FeatureExtraction_Association_PointToPointAssociation

<!-- PROJECT LOGO -->
<br />
<p align="center">
  <!-- <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a> -->

  <h3 align="center">FeatureExtraction_Association_PointToPointAssociation</h3>

  <p align="center">
    MATLAB code implementation of a series of functions that associate spatial
data. Specifically, functions are provided to determine matches between data
sets of (X,Y) points, store and compare groups of associated points (patch
objects), and determine intersections between patch objects and circular arcs
(useful for collision detection).
    <br />
    <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/tree/main/Documents">View Demo</a>
    ·
    <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/issues">Report Bug</a>
    ·
    <a href="https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/issues">Request Feature</a>
  </p>
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
data. Specifically, functions are provided to determine matches between data sets of (X,Y) points, store and compare groups of associated points (patch objects), and determine intersections between patch objects and circular arcs (useful for collision detection).

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
The majority of the code for the point and patch association functionalities are implemented in this directory. All functions as well as test scripts are provided.

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
<ul>
datasets: a structure array containing subfields of X and Y 
coordinates in the following form:

    datasets{i_set}.X   
    datasets{i_set}.Y
</ul>
Note that i_set denotes a data set structure. Each set can be
plotted separately and all the dataset can also be plotted in the same figure.  

**OUTPUTS:**

h: a handle to the resulting figure

**Examples:**

The examples for this function are present in **script_test_fcn_Points_plotSetsXY.m**

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
%% Setup

% Call the function to fill in an array of "dataset" type
datasets_array = fcn_Points_fillPointSampleSets;


%% Example 1
% We can save one of these as a single "dataset"

single_path = {datasets_array{1}};

fcn_Points_plotSetsXY(single_path,1);
text(single_path{1}(1,1),single_path{1}(1,2),'Start');
```
<pre align="center">
  <img src=".\Images\fcn_Points_fillPointSampleSets_Ex1.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Plot for the 1st dataset in fcn_Points_fillPointSampleSets.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

```MATLAB

%% Example 2
% illustrating the different XY points on different plots

for i_set = 1:length(datasets_array)
    fcn_Points_plotSetsXY(datasets_array(i_set),(i_set+1));
    text(datasets_array{i_set}(1,1),datasets_array{i_set}(1,2),'Start');
end
```
<pre align="center">
  <img src=".\Images\fcn_Points_fillPointSampleSets_Ex2_1.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Plot for the 1st dataset in fcn_Points_fillPointSampleSets.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<pre align="center">
  <img src=".\Images\fcn_Points_fillPointSampleSets_Ex2_2.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Plot for the 2nd dataset in fcn_Points_fillPointSampleSets.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>

<pre align="center">
  <img src=".\Images\fcn_Points_fillPointSampleSets_Ex2_3.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Plot for the 3rd dataset in fcn_Points_fillPointSampleSets.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
</pre>
```MATLAB
%% Example 3
% All plots on the same graph

fcn_Points_plotSetsXY(datasets_array,5);

```

<pre align="center">
  <img src=".\Images\fcn_Points_fillPointSampleSets_Ex3.jpg" alt="fcn_Points_fillPointSampleSets picture" width="400" height="300">
  <figcaption>Plot for the all the datasets in fcn_Points_fillPointSampleSets.</figcaption>
  <!--font size="-2">Photo by <a href="https://unsplash.com/ko/@samuelchenard?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Samuel Chenard</a> on <a href="https://unsplash.com/photos/Bdc8uzY9EPw?utm_source=unsplash&utm_medium=referral&utm_content=creditCopyText">Unsplash</a></font -->
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

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### fcn_Points_adjustPointSetStatistics

Adjusts the XY positions of all datasets existing in a data structure by
the provided values

FORMAT: 

```MATLAB
      h = fcn_Points_adjustPointSetStatistics(datasets,bias,noiseMean,noiseVariance)
```
**INPUTS:**

     datasets: a structure array containing subfields of X and Y 
     coordinates in the following form:
          datasets{i_set}.X
          datasets{i_set}.Y
     Note that i_set addresses a particular data set structure. Each set
     will be modified separately.
     bias: a matrix that defines constant X and Y offsets for each data 
           set, hence an N x 2 matrix where N is the number of data sets
     noiseMean: a matrix that defines the X and Y means for the noise
           to be added to each of the each data sets, hence N x 2
     noiseVariance: a matrix that defines the X and Y variance for the
           noise to be added to each of the each data sets, hence N x 2

**OUTPUTS:**

     datasets: the input datasets will be modified and returned

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
          datasets{i_set}.X
          datasets{i_set}.Y
     Note that i_set addresses a particular data set structure. Each set
     will be modified separately.
     R: radius of max noise

                  OPTIONAL INPUTS:

     fig_num: figure number

**OUTPUTS:**

     datasets_out: the input datasets will be modified and returned


<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### fcn_Points_calcPairStatistics

Calculates the statistics of the errors between two sets of paired XY
data points and an overall X,Y map shift to match the centroids of the
data sets.


FORMAT:

      [errRMS,errVar,meanShift] = fcn_Points_calcPairStatistics(pairedXYdata)

INPUTS:

     pairedXYdata: an N x 4 matrix with [X1 Y1 X2 Y2] data in each row.

OUTPUTS:

     errRMS: the RMS distances associated with the point-to-point match errors
     errVar: the variance in the point-to-point match error distances
     meanShift: the [X Y] shift that would cause the centroids of the
                two data sets to match

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

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### script_AssociatePointsIncremental

This script simulates a vehicle following a trajectory by detecting the nearest points usings KNN algoorithm.


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

