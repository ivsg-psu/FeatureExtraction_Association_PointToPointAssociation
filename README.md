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
      <a href="#structure">Structure</a>
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
          <li><a href="#fcn_points_checkinputstofunctions">fcn_Points_checkInputsToFunctions</a></li>
          <li><a href="#fcn_points_fillpointsamplesets">fcn_Points_fillPointSampleSets</a></li>
          <li><a href="#fcn_points_fillpointsetviauserinputs">fcn_Points_fillPointSetViaUserInputs</a></li>
          <li><a href="#fcn_points_plotsetsxy">fcn_Points_plotSetsXY</a></li>
          <li><a href="#fcn_points_pairxydata">fcn_Points_pairXYdata</a></li>
          <li><a href="#fcn_points_calcpairstatistics">fcn_Points_calcPairStatistics</a></li>
          <li><a href="#fcn_points_adjustpointsetstatistics">fcn_Points_adjustPointSetStatistics</a></li>
        </ul>
        </li>
        <li><a href="#patch-object-creation-and-manipulation-functions">Patch Object Creation and Manipulation Functions</a></li>
        <ul>
          <li><a href="#fcn_patch_fillsamplepatches">fcn_Patch_fillSamplePatches</a></li>
          <li><a href="#fcn_patch_fillpatcharrayviauserinputs">fcn_Patch_fillPatchArrayViaUserInputs</a></li>
          <li><a href="#fcn_patch_plotpatch">fcn_Patch_plotPatch</a></li>
          <li><a href="#fcn_patch_insertpoints">fcn_Patch_insertPoints</a></li>
          <li><a href="#fcn_patch_determineaabb">fcn_Patch_determineAABB</a></li>
          <li><a href="#fcn_patch_inferprimitive">fcn_Patch_inferPrimitive</a></li>
          <li><a href="#fcn_patch_checkcollisions">fcn_Patch_checkCollisions</a></li>
        </ul>
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
<!-- FUNCTION DEFINITIONS -->
## Functions
The majority of the code for the point and patch association functionalities are implemented in this directory. All functions as well as test scripts are provided.

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

### Point-Set Association Functions
<ul>
	<li>fcn_Points_checkInputsToFunctions: TEMPLATE function for checking arguments to functions, such as point sets, etc. to make sure the formatting and sizes are correct</li>
	<li>fcn_Points_fillPointSampleSets: a function to load some sample data sets to use for testing the other functions</li>
	<li>fcn_Points_fillPointSetViaUserInputs: a function that allows a user to create (X,Y) point sets by clicking in a figure with the mouse</li>
	<li>fcn_Points_plotSetsXY: a function that plots (X,Y) point sets with various options</li>
	<li>fcn_Points_pairXYdata: a function that associates the mutually closest points in two different point sets and returns the pairs as well as points which don't have an obvious mutual pair, in both directions</li>
	<li>fcn_Points_calcPairStatistics: a function that calculates the statistics for paired sets of points, returning RMS deviation, variance in point locations, and the offset between the centroids of the two point sets (a measurement of the systematic "shift" between two point sets)</li>
	<li>fcn_Points_adjustPointSetStatistics: a function to add 2D Gaussian noise and/or bias to a point set (e.g. to simulate sensor noise or bias) </li>
</ul>

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>

***

#### fcn_Points_checkInputsToFunctions

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### fcn_Points_fillPointSampleSets

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### fcn_Points_fillPointSetViaUserInputs

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### fcn_Points_plotSetsXY

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### fcn_Points_pairXYdata

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### fcn_Points_calcPairStatistics

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### fcn_Points_adjustPointSetStatistics

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

### Patch Object Creation and Manipulation Functions
<ul>
	<li>fcn_Patch_fillSamplePatches: a function to load a few sample patch objects of different sizes, shapes, and colors for testing the other patch object functions</li>
	<li>fcn_Patch_fillPatchArrayViaUserInputs: a function that allows a user to create patch objects by choosing patch colors and then using a mouse to click in a figure window to define the vertices of the patch object</li>
	<li>fcn_Patch_plotPatch: a function that plots patch a patch object or an array of patch objects, optionally choosing particular patch objects from the array and/or plotting into a particular figure</li>
	<li>fcn_Patch_insertPoints: a function that allows the user to insert one or more (X,Y) points into a patch object to add vertices to the patch object</li>
	<li>fcn_Patch_determineAABB: a function to determine the (X,Y) extremes of the patch object and store them in the patch object attributes</li>
	<li>fcn_Patch_inferPrimitive: a function to test the fit of circular and rectangular shape primitives to the vertices of the patch object and store the best fit to the patch object attributes (or reject both fits and label the object as irregular)</li>
	<li>fcn_Patch_checkCollisions: a function to test an array of patch objects to determine impact point and time or, for non-collisions, the closest point and time of closest approach with a rectangular object traveling in a circular trajectory of a given radius, center point, and initial heading</li>
</ul>

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>

***

#### fcn_Patch_fillSamplePatches

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### fcn_Patch_fillPatchArrayViaUserInputs

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### fcn_Patch_plotPatch

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### fcn_Patch_insertPoints

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### fcn_Patch_determineAABB

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### fcn_Patch_inferPrimitive

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***

#### fcn_Patch_checkCollisions

<a href="#featureextraction_association_pointtopointassociation">Back to top</a>
***
Each of the functions has an associated test script, using the convention
	```sh
	script_test_fcn_fcnname
	```
where fcnname is the function name starting with "Patch" or "Point" as listed above.


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



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

Your Name - [@twitter_handle](https://twitter.com/twitter_handle) - email

Project Link: [https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation](https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation)



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

