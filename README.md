# Beaufort-Sea-Lagoons
Supplementary Materials for Harris et al. (2018)
Do high Arctic coastal food webs rely on a terrestrial carbon subsidy?
Food Webs

In RStudio, open Beaufort_lagoons_simmr.R to run analysis. Be sure all files from this repository are in your working directory.

JAGS must be installed on your computer to run some of these scripts. Install it here: https://sourceforge.net/projects/mcmc-jags/files/

Some of the analysis scripts have been borrowed from other sources:

  McTigue & Dunton (2017) --> see https://github.com/nathanmct/Chukchi-Sea-trophodynamics

  From Andrew Parnell's GitHub for simmr (https://github.com/andrewcparnell/simmr/tree/master/R). The functionality is the same, but I have altered some of the ggplot code to change the output visualization.

Much of it has been modified for our own purposes, but we acknowledge and thank these sources for their open-source code. The script outlines the analysis and visualization of the stable isotope analysis in the paper.

File descriptions:
BeaufortSeaLagoons_isotope_data.csv is a file containing all stable isotope data for the study.

source_averages.csv contains the end-members for the invertebrates and fish.
source.averages.mammals.csv contains the end-members for mammals.

corrections.csv contains the trophic enrichment factors for each end-member for Su/FF, Ss/De, Ep/Om.
corrections.higherTL.csv contains the trophic enrichment factors for each end-member for Fish.
corrections.mammTL.csv contains the trophic enrichment factors for each end-member for Mam/Carn.

The following are files used to place rectangles and labels in Fig. 4
MPB_polygons.csv
MPB_Label.csv
Shelf_polygons.csv
POM_Label.csv
Terrestrial_Label.csv
Terrestrial_polygons.csv
