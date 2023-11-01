# pupilArousal
Pupillometry project experimental and analysis scripts for MATLAB.

## Authors ########################################################
June Hee Kim &lt;<junehee.kim@nih.gov>&gt; Zvi N. Roth &lt;<zviroth@gmail.com>&gt; and Elisha P. Merriam &lt;<elisha.merriam@nih.gov>&gt; 

## Usage ##########################################################
Experimental and core analysis scripts used for Kim J, Yin C, Merriam EP, Roth ZN (2023) Pupil Size Is Sensitive to Low-Level Stimulus Features, Independent of Arousal-Related Modulation. ENeuro 

Scripts can be roughly broken down as follows. Scripts are listed from innermost function to outermost (i.e., the first function is called by the second function, the second is called by the third, etc)

1. **Experimental Code**
   rwdRapid.m: This code measures arousal and stimulus effects on pupil responses by executing a rapid serial visual presentation (RSVP) task. During this task, numerical digits swiftly succeed one another at the center of the screen. Each of the 17 trials in a run lasts for 15 seconds. While each run corresponds to either a high or low reward level, the initial run serves as a practice session. In this practice run, the stimulus onset asynchrony (SOA) is set using a one-up-two-down staircase method, which then informs the SOA for subsequent runs. Note: setup of MGL is required to run this code. Please refer to: https://gru.stanford.edu/doku.php/mgl/overview
   
3. **Data Extraction & Pre-processing**

   Utilities: myGetTaskEyeTraces.m, myMglEyelinkEDFRead.m,  extractPupilData.m (used to read eye data (.edf) and corresponding stim files (.mat) for each run. Blink interpolation is handled with data extraction).

   paramorg.m: Filters extracted dataset to check for exclusion criteria, variables, and recorded data.
   
4. **Analysis**

   bootstrapPA.m: Validates the linear regression analysis used to calculate response amplitudes.
   
   pupilArousal.m: This is an analysis script designed to evaluate task-related and stimulus-evoked pupil responses. The script undertakes several specific operations, including:

   * Z-scoring pupil sizes within individual runs,
   * Implementing a low-pass Butterworth filter (3rd order) on the pupil data,
   * Constructing stimulus and null response templates,
   * Estimating trial-wise response amplitudes using linear regression.

   Furthermore, the accuracy of our regression analysis for determining response amplitude has been validated through a bootstrapping process.
   
## Example Data ###################################################
We include real data from two subjects in our experiment to show how the input and output data are structured for the analysis script. PupilData consists of deanonymized data originally extracted using extractPupilData.m with blinks interpolated in the process.

PupilData{N, 1}{run} -- extracted eye data(position, diameter, reaction time, etc.), visual stimuli values presented during tasks, and associated labels for each stimuli value. 

PupilData{N, 2}{run} -- correctness and sample rate

PupilData{N, 3}{run} -- task variable information for each trial (spatial frequency, contrast, and null trials) 

PupilData{N, 4}{run} -- reward level

## Dependencies ###################################################
Statistics and Machine Learning Toolbox, MGL

## Citing #########################################################

To cite, please reference the following: June Hee Kim, Christine Yin, Elisha P. Merriam, Zvi N. Roth (2023) Pupil Size Is Sensitive to Low-Level Stimulus Features, Independent of Arousal-Related Modulation. eneuro 10:ENEURO.0005-23.2023.

## References #####################################################

* de Gee JW, Colizoli O, Kloosterman NA, Knapen T, Nieuwenhuis S,Donner TH (2017) Dynamic modulation of decision biases by brainstem arousal systems. Elife 6:e23232.
* Gardner JL, Merriam EP, Schluppeck D, Larsson J (2018) MGL: visual psychophysics stimuli and experimental design package. DOI (10.5281/zenodo.129949)
* Levitt H (1971) Transformed Up‐Down Methods in Psychoacoustics. J Acoust Soc Am 49:467–477.
* Roth ZN, Ryoo M, Merriam EP (2020) Task-related activity in human visual cortex. PLoS Biol 18:e3000921.

