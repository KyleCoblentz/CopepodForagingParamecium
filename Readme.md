# Readme

This repository contains the data and code for Coblentz, K.E., L. Yang, A. Dalal, M. Incarnato, D. Thilakarathne, C. Shaw, R. Wilson,
F. Biagioli, K.L. Montooth, and J.P. DeLong. 2024. Heritable intraspecific variation among prey in size and movement interact to shape
predation risk and potential natural selection. Function Ecology.

### Data

The repository contains the following .csv files containing the raw data from the experiment:

CopepodForagingLengths_Clean.csv -- This file contains three columns: *Copepod.ID* which contains the label for each copepod, *Length* which gives
each copepods length in mm, and *Width* which gives the width of each copepod in mm.

Paramecium_Phenotype_Data.csv -- This file contains the morphological and movement data for each cell maeasured across the 126 outcrossed 
*Paramecium* lines in the experiment. The final column contains a variable *id* that has a string for each row, and the component of the 
string consisting of "G_[Number]" gives the name of the outcrossed line that the cell belongs to. For explanations of the 
morphological and movement variables, see the supplementary material of the published manuscript.

Foraging trial data.csv -- This file contains the data from foraging trials of copepods consuming *Parmamecium*. The file contains information
on the date on which the foraging trial was performed (*Date*), the ID of the of individual copepod used in the trial (*Copepod ID*), the ID number of 
the outcrossed line used in the foraging trial (*Prey type*, note that this corresponds to the number in the G_[Number] string in the id column of the 
phenotype data), the time at which the trial was started and ended (*Time strarted* and *Time ended*), the prey offered at the beginning of the trial and 
the prey remaining at the end of the trial (*Prey offered* and *Prey remaining*), a column for notes (*Notes*), and a column designating the 
copepod that died during the experiment (*Dead*).

## Code

The repository contains the two following .R files that were used to analyze the data:

QuantitativeGeneticsAnalysis_Clean.R -- This R file contains the code to quantify the heritability of the *Paramecium* traits from
the *Paramecium_Phenotype_Data.csv* file.

CopepodForagingAnalysis_Clean.R -- This R file contains the code to perform a princpal components analysis on the *Paramecium* phenotype 
data, perform the regression analysis on the proportion of *Paramecium* consumed during the foraging trials, and generates all of the figures 
within the the manuscript and supplementary material.

