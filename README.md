
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Thanks again! Now go create something AMAZING! :D
***
***
***
*** To avoid retyping too much info. Do a search and replace for the following:
*** github_username, repo_name, twitter_handle, email, project_title, project_description
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]



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
    <li><a href="structure">Repo Structure</a>
	    <ul>
	    <li><a href="#directories">Top-Level Directories</li>
	    <li><a href="#functions">Functions</li>
	    </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

<!--[![Product Name Screen Shot][product-screenshot]](https://example.com)-->

MATLAB code implementation of a series of functions that associate spatial
data. Specifically, functions are provided to determine matches between data sets of (X,Y) points, store and compare groups of associated points (patch objects), and determine intersections between patch objects and circular arcs (useful for collision detection).



<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation.git
   ```

<!-- STRUCTURE OF THE REPO -->
## Structure
### Directories
The following are the top level directories within the repository:
<ul>
	<li>Documents: Descriptions of the functionality and usage of the various MATLAB functions and scripts in the repository.</li>
	<li>Functions: The majority of the code for the point and patch association functionalities are implemented in this directory. All functions as well as test scripts are provided.</li>
	<li>Utilities: Dependencies that are utilized but not implemented in this repository are placed in the Utilities directory. These can be single files but are most often other cloned repositories.</li>
</ul>

<!-- FUNCTION DEFINITIONS -->
### Functions
**Point-Set Association Functions**
<ul>
	<li>fcn_Points_checkInputsToFunctions: TEMPLATE function for checking arguments to functions, such as point sets, etc. to make sure the formatting and sizes are correct</li>
	<li>fcn_Points_fillPointSampleSets: a function to load some sample data sets to use for testing the other functions</li>
	<li>fcn_Points_fillPointSetViaUserInputs: a function that allows a user to create (X,Y) point sets by clicking in a figure with the mouse</li>
	<li>fcn_Points_plotSetsXY: a function that plots (X,Y) point sets with various options</li>
	<li>fcn_Points_pairXYdata: a function that associates the mutually closest points in two different point sets and returns the pairs as well as points which don't have an obvious mutual pair, in both directions</li>
	<li>fcn_Points_calcPairStatistics: a function that calculates the statistics for paired sets of points, returning RMS deviation, variance in point locations, and the offset between the centroids of the two point sets (a measurement of the systematic "shift" between two point sets)</li>
	<li>fcn_Points_adjustPointSetStatistics: a function to add 2D Gaussian noise and/or bias to a point set (e.g. to simulate sensor noise or bias) </li>
</ul>

**Patch Object Creation/Manipulation Functions**
<ul>
	<li>fcn_Patch_fillSamplePatches: a function to load a few sample patch objects of different sizes, shapes, and colors for testing the other patch object functions</li>
	<li>fcn_Patch_fillPatchArrayViaUserInputs: a function that allows a user to create patch objects by choosing patch colors and then using a mouse to click in a figure window to define the vertices of the patch object</li>
	<li>fcn_Patch_plotPatch: a function that plots patch a patch object or an array of patch objects, optionally choosing particular patch objects from the array and/or plotting into a particular figure</li>
	<li>fcn_Patch_insertPoints: a function that allows the user to insert one or more (X,Y) points into a patch object to add vertices to the patch object</li>
	<li>fcn_Patch_determineAABB: a function to determine the (X,Y) extremes of the patch object and store them in the patch object attributes</li>
	<li>fcn_Patch_inferPrimitive: a function to test the fit of circular and rectangular shape primitives to the vertices of the patch object and store the best fit to the patch object attributes (or reject both fits and label the object as irregular)</li>
	<li>fcn_Patch_checkCollisions: a function to test an array of patch objects to determine impact point and time or, for non-collisions, the closest point and time of closest approach with a rectangular object traveling in a circular trajectory of a given radius, center point, and initial heading</li>
</ul>
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
_For more examples, please refer to the [Documentation] https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation/tree/main/Documents)_



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

Your Name - [@twitter_handle](https://twitter.com/twitter_handle) - email

Project Link: [https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation](https://github.com/ivsg-psu/FeatureExtraction_Association_PointToPointAssociation)




<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/ivsg-psu/repo.svg?style=for-the-badge
[contributors-url]: https://github.com/ivsg-psu/repo/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/ivsg-psu/repo.svg?style=for-the-badge
[forks-url]: https://github.com/ivsg-psu/repo/network/members
[stars-shield]: https://img.shields.io/github/stars/ivsg-psu/repo.svg?style=for-the-badge
[stars-url]: https://github.com/ivsg-psu/repo/stargazers
[issues-shield]: https://img.shields.io/github/issues/ivsg-psu/repo.svg?style=for-the-badge
[issues-url]: https://github.com/ivsg-psu/repo/issues
[license-shield]: https://img.shields.io/github/license/ivsg-psu/repo.svg?style=for-the-badge
[license-url]: https://github.com/ivsg-psu/repo/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/ivsg-psu

