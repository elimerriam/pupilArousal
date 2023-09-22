# pupilArousal
Pupillometry project analysis scripts for MATLAB.

## Authors ########################################################
June Hee Kim &lt;<junehee.kim@nih.gov>&gt; Zvi N. Roth &lt;<zviroth@gmail.com>&gt; and Elisha P. Merriam &lt;<elisha.merriam@nih.gov>&gt; 

## Usage ##########################################################
Core analysis scripts used for Kim J, Yin C, Merriam EP, Roth ZN (2023) Pupil Size Is Sensitive to Low-Level Stimulus Features, Independent of Arousal-Related Modulation. ENeuro 

Scripts can be roughly broken down as follows. Scripts are listed from innermost function to outermost (i.e., the first function is called by second function, second is called by third etc)

1. ** Data Extraction & Pre-processing **

   Utilities: myGetTaskEyeTraces.m, myMglEyelinkEDFRead.m

   extractPupilData.m: Used to read eye data (.edf) and corresponding stim files (.mat) for each run. Blink interpolation is handled with data extraction.

   paramorg.m: Filters extracted dataset to check for variables and recorded data.
   
3. **Analysis**

   bootstrapPA.m: Validates the linear regression analysis used to calculate response amplitudes.
   
   pupilArousal.m: Analysis script for task-related and stimulus-evoked pupil responses. Specifically does the following: z-score pupil sizes within each run, low pass Butterworth (3rd order) filtering of pupil data, calculating stim & null response template, and estimating trial wise response amplitude through linear regression. Our regression analysis for estimating response amplitude has been evaluated using bootstrapping procedure.
   
## Example Data ###################################################
We include real data from two subjects in our experiment, so you can see how the input and output data are structured for the analysis script. PupilData consists of deanonymize data originally extracted using extractPupilData.m with blinks interpolated in the process. 

PupilData{N, 1} -- extracted eye(position, diameter, saccade, etc..) by viable run

PupilData{N, 2} -- stimulus, screen, and task info by viable run 

PupilData{N, 3} -- task variable information by viable run 

PupilData{N, 4} -- reward info by viable run

## Dependencies ###################################################
Statistics and Machine Learning Toolbox, MGL

## Citing #########################################################

To cite, please reference the following:

## References #####################################################

* de Gee JW, Colizoli O, Kloosterman NA, Knapen T, Nieuwenhuis S,Donner TH (2017) Dynamic modulation of decision biases by brainstem arousal systems. Elife 6:e23232.
* Gardner JL, Merriam EP, Schluppeck D, Larsson J (2018) MGL: visual psychophysics stimuli and experimental design package. DOI (10.5281/zenodo.129949)
* Roth ZN, Ryoo M, Merriam EP (2020) Task-related activity in human visual cortex. PLoS Biol 18:e3000921.

