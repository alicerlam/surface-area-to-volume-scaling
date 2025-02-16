# Overview
A MATLAB workflow for analyzing cell line profiles in light microscopy images in order to quanitfy plasma membrane volume relative to cell size. This collection of scripts allows users to select one slice of interest of a fluorescent Z-stacked timelapse movie for each cell identified in the original image, before automatically calculating the radius of the cell, the background intensity inside and outside of the cell, and the maximum intensity/full-width-half-maximum intensity/area of the cell membrane.

# Installation:
1. Download the repo using the green "Clone or download" button on this page.
2. Alternatively: Make sure git is installed and in the terminal enter `git clone https://git@github.com/alicerlam/surface-area-to-volume-scaling.git`
3. Depending on your file format, you may also need to download bfmatlab (follow the installation instructions here: https://docs.openmicroscopy.org/bio-formats/6.3.1/users/matlab/index.html). bfmatlab was used in our paper to handle .dv files.

# Usage
Open select_images.m in MATLAB2023a. 
1. Under the section labeled "Run me ONCE per experiment/folder!!", set your paths to the directories containing bfmatlab (if needed) and your images under the first comment.
2. Set variable "dZ" to reflect the number of slices within your Z stack, and "dT" to reflect the number of timepoints in your timelapse data.
3. Set "output_file" to a location for the output spreadsheet.
4. To account for cell movement during the timelapse, there is an optional late_timepoint variable such that you can reselect the Z-slice of interest at one later timepoint and downstream analysis will utilise the later selected Z-slice. Both Z-slice indices are saved. If the cells are stationary, it can be set to 1.
5. Under "Setup", find the find_files function and change the argument to match your image names. find_files is capable of taking wildcards in order to select all images with a given substring in the name.
6. Under "Parameters", set the start and end image number to match the number of images in your image directory.
7. The code will attempt to automatically identify the cells using a Hough transform. You can lower the "sensitivity" variable if it is identifying things that are not cells, and vice versa. If the automatic segmentation does not work well, set variable "usemanual" to 1 or answer 0 when prompted "Did I find the right number of circles? Enter 1 if yes, 0 if no.". At that point, cells can be manually selected using the CROIEditor GUI that pops up.
8. For each cutout of a cell at a timepoint, enter a list of Z-slice indices corresponding to the slices of interest.
9. If late_timepoint is set to a timepoint other than 1, steps 7-8 will repeat at the timepoint set.
10. The output file contains one sheet containing the cell radius, the background intensity inside and outside of the cell, and the maximum intensity/full-width-half-maximum intensity/area of the cell membrane for each timepoint in the folder. The last sheet contains data for which Z-stacks were selected for analysis at the initial and later timepoints, as well as the image names.

# What does this code do?
1. Cells are roughly identified either by a Hough transform (automatic) or by manual selection with CROIEditor. ![alt_text](https://i.imgur.com/Vsa0hdf.jpeg)
2. For every curated cell, there is an image that has all cells other than the cell of interest blacked out to prevent line profile collisions.![alt text](https://i.imgur.com/1bPc4aq.jpeg)
3. 360 line profiles are taken from the center of the cell moving radially outward. The cell membrane is identified by a peak in pixel intensity, shown in red below. ![alt text](https://i.imgur.com/MTHVZ4R.jpeg)
4. Line profiles are automatically aligned to the peak in pixel intensity and assumes an approximately circular cell. Plotting the line profiles can identify which cells should be excluded from further analysis due to poor peak alignment (e.g. Cell 4). ![alt text](https://i.imgur.com/BCoszOy.png)
5. Output statistics can be further analysed to investigate changes in membrane volume ('AUC') as a function of cell radius, time, or other perturbation. An example analysis of membrane volume by radius (comparing membrane volume in small and large cells) and how the relative change evolves over time is provided below.
   ![alt text](https://i.imgur.com/MfEfECq.png) ![alt text](https://i.imgur.com/V80WTgi.png)
   
# Cite as:
Wu W, Lam AR, Suarez K, Smith GN, Duquette SM, Yu J, Mankus D, Bisher M, Lytton-Jean A, Manalis SR, Miettinen TP. Constant surface area-to-volume ratio during cell growth as a design principle in mammalian cells. bioRxiv [Preprint]. 2024 Jul 18:2024.07.02.601447. doi: 10.1101/2024.07.02.601447. PMID: 39005340; PMCID: PMC11244959.\
\
**CROIEditor.m:** Jonas Reber (2025). Multi ROI/Mask Editor Class (https://www.mathworks.com/matlabcentral/fileexchange/31388-multi-roi-mask-editor-class), MATLAB Central File Exchange. Retrieved February 16, 2025.
