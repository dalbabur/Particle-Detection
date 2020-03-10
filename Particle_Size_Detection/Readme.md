# Particle Size Detection Within Vortex

The aim of this subfolder is to determine the size distribution of particles that are trapped within the vortex. The requirement is that you have focused videos, where the particle is edge is clearly defined, and the file is in a tiff format (multiple or single page is fine). **IF** multipage, then you should know exactly which frame you want to edit.

**WORKFLOW FOR ANALYIS**
1. Take videos using high speed camera, with frame rate > particle speed so that the images are still
2. Convert segment of interest from .cine to .tif (or .tiff, but you need to make additional code changes)
3. Use FIJI/ImageJ to determine which frames have the most particles in focus with distinct edges
4. Open and run program titled 'bead_size_measurements_6um_beads'
5. Follow prompts, output file will be an excel spreadsheet

When opening the Excel spreadsheet, the order of the columns will be the following (Left to right)
* Particle Detection number
* Center location (X) in MATLAB image coordinates
* Center location (Y) in MATLAB image coordinates
* **Radius** of particle **(IN PIXELS)**
* Found by the algorithm, identified with human intervention, or not a particle
  * 1 - found by the algorithm
  * 3 - Human detection
  * 0 - Determined not to be a new particle

The program uses 2 filters - one is a Sobel filter, and one is a tophat filter. Both work for edge detection, and you can modify the parameters to see if there's a better combination of parameters for your focusing parameters. When you download, make sure you add the path of the GSobel Filter to the MATLAB Path

To validate the code, we use the inlet conditions and the ground truth data from the Scepter (handheld cell counter in lab) to ensure that our code correctly measures the particle size with respect to the Scepter.

As you go through the data, you should know the following:
* The box that it asks you to draw, should be drawn only around the particles you wish to analyze in that particular frame
  * If particles are off focus in that frame but are focused in other frames, I usually won't draw the box around the particles that are    out of focus, and find other frames where they are focused
* The pause after it has identified the particles is to move the particle window so that you can look at the image while evaluating your   options

If you would like to plot the output of the scepter data, use the following: Scepter_freq_2_data.m. This code takes in the second column (user has to chose the file) of the scepter data, and calculates the distribution. It then gates the data by the smallest expected cell size.
