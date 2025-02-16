# Overview
A collection of scripts that allows users to select one slice of interest of a fluorescent Z-stacked timelapse movie for each cell identified in the original image, before automatically calculating the radius of the cell, the background intensity inside and outside of the cell, and the maximum intensity/full-width-half-maximum intensity/area of the cell membrane.

# Installation:
1. Download the repo using the green "Clone or download" button on this page.
2. Alternatively: Make sure git is installed and in the terminal enter `git clone https://git@github.com/alicerlam/surface-area-to-volume-scaling.git`
3. Depending on your file format, you may also need to download bfmatlab (follow the installation instructions here: https://docs.openmicroscopy.org/bio-formats/6.3.1/users/matlab/index.html). bfmatlab was used in our paper to handle .dv files.

# Usage
Open select_images.m in MATLAB2023a. Under the section labeled "Run me ONCE per experiment/folder!!", set your paths to the directories containing bfmatlab (if needed) and your images under the first comment. Set variable "dZ" to reflect the number of slices within your Z stack, and "dT" to reflect the number of timepoints in your timelapse data.

# Cite as:
Wu W, Lam AR, Suarez K, Smith GN, Duquette SM, Yu J, Mankus D, Bisher M, Lytton-Jean A, Manalis SR, Miettinen TP. Constant surface area-to-volume ratio during cell growth as a design principle in mammalian cells. bioRxiv [Preprint]. 2024 Jul 18:2024.07.02.601447. doi: 10.1101/2024.07.02.601447. PMID: 39005340; PMCID: PMC11244959.\
\
**CROIEditor.m:** Jonas Reber (2025). Multi ROI/Mask Editor Class (https://www.mathworks.com/matlabcentral/fileexchange/31388-multi-roi-mask-editor-class), MATLAB Central File Exchange. Retrieved February 16, 2025.
