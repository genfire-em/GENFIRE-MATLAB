#_Welcome to GENFIRE_

Welcome to the MATLAB implementation of GENFIRE. 


GENFIRE is a robust, Fourier-based reconstruction algorithm that is
capable of using a limited set of input projections to generate a 3D reconstruction
while also partially retrieving missing projection information. It does this by iterating 
between real and reciprocal space and applying simple constraints in each to find an optimal solution that  
correlates with the input projections while simultaneously obeying real space conditions
such as positivity and support.The result is a more consistent and faithful reconstruction with superior contrast and, in many cases, resolution when
compared with more traditional 3D reconstruction algorithms such as Filtered Back-Projection (FBP) or the Algebraic Reconstruction Technique (ART).  

## Useage
Reconstructions are run either with `GENFIRE_Main.m` for a cubic array and 3 Euler angles or with `GENFIRE_Main_Tomo.m` for a rectangular array with a single tilt-axis. Simply edit the parameters in the appropriate file and run the script to execute the reconstruction. The parameters that may be adjusted are

* filename_Projections (str) - filename containing projection images as an N x N x num\_projections array

* filename_Angles (str) - filename containing the angles as either a num\_projections x 3 array of Euler angle triples (phi, theta, psi) or a num\_projections x 1 array indicating a single-axis tilt series
* filename_Support (str) - filename containing an NxNxN binary support
* numIterations (int) - number of GENFIRE iterations to run
* pixelSize (double) - size of a single pixel, used to display the Fourier Shell Correlation (FSC) with proper units if desired. Setting `pixelSize` to 0.5 will 
* oversamplingRatio (int) - controls the amount of zero padding. Formally, the oversampling ratio in a given direction is the total array size (after padding) divided by the size of the input projection.
* griddingMethod (int) - controls the gridding method. Choose from the following:
	1. FFT - Faster, less accurate
	2. DFT - Slower, more accurate
* constraintEnforcementMode (int) - choose a method of enforcing the Fourier constraint from:
	1. Resolution extension and suppression
	2. Resolution extension only
	3. Enforce all datapoints always (No resolution extension or suppression)
* interpolationCutoffDistance (double) - maximum tolerable radial distance to consider measured datapoints for gridding to a given Fourier voxel

* ComputeFourierShellCorrelation (bool) - the data will be divided into two halves, independently reconstructed, and the FSC will be compared between the two halves. A final reconstruction will then be computed from all of the data 
* numBins (int) - number of bins to use in the calculation of FSC (if applicable)  
* percentValuesForRfree - percentage of data to withold within each shell for calculation of Rfree
* numBinsRfree - number of bins to divide Fourier space in for Rfree
* doCTFcorrection (bool) *(experimental)* - turn on CTF correction
* CTFThrowOutThreshhold (double) *(experimental)* - Fourier values in regions where the corresponding CTF is lower than this value will not be gridded
