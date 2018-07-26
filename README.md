# Matlab based GMMUBM and Identity-Vector Classification for EEGs:

The files in this repository are provided to enable interested readers to perform the experiments outlined in the paper published in The Journal of Neuroscience Methods.

The only GMM-UBM and I-Vector tool that existed at the start of my research came from from the [MSR Identitiy Toolbox] (http://research.microsoft.com/en-us/downloads/2476c44a-1f63-4fe0-b805-8c2de395bb2c/). However, after experimenting with the MSR Toolbox I found it lacking. I then built this framework to better suit my needs in terms of ease of use and compliance with newer versions of Matlab.

All of the work is created by me, <christian.radcliffe.ward@gmail.com>, during my tenure at Temple University working under the umbrella of the [NEDC](https://www.isip.piconepress.com/projects/nedc/) under the guidance of [Dr. Iyad Obeid](https://engineering.temple.edu/faculty/iyad-obeid).

### Caution

These files were written to be run on a computing cluster operating **qsub**. They take advantage of parfor loops and do not use any gpuarrays. If run on a local machine it may crash matlab or the local machine. Be careful.

# Functionality:
1. Process raw data (EDF) into features
2. Organize features into epochs
3. Evaluate epochs through Mahalanobis, GMMUBM, and I-Vector algorithms
4. Organize evaluations into results

## Process raw data

To produce the cepstrum (HTK), power spectral density (PSD), and spectral coherence (COH) features the following functions are used to modify the EDF files.

The PSD and COH features were based on the work of La Rocca et al's 2014 paper: _Human Brain Distinctiveness Based on EEG Spectral Coherence Connectivity_ 

 ```
 laRoccaFeatureBuilder(file_list, feature_type, final_folder, n_epochs, workers)
 ```
 
 Given the a **file_list** of EDF files and a desired **feature_type** (0 for PSD, 1 for COH) it will partition the data of each channel into **n_epochs**. They will be stored in **final_folder** (full path name) and be carried out by the number of specified **workers**. While written for the PhysioNet EEGMMIDB, the functionality should extend to any set of EDF files.
 
 To produce the CEP features, an external toolbox is required: [VOICEBOX](http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html). This kit, written by Mike Brookes, provides tools to convert the EDF signals to cepstrum htk files and then read those .HTK files. This is original a tool for speech processing, but why re-invent what already  exists.
 
 To produce the HTK features the same function is used, but **feature_type** is set to 2 for CEP. This approach requires the incoming list to be a list of .HTK files from the already processed .EDF files. Producing the HTK features, 1 second windows with 1/10 second frames using a TCP montage resulting in 22 channels, is outlined in a [paper](https://arxiv.org/ftp/arxiv/papers/1801/1801.02477.pdf) by Harati et al.
 
 ```
 laRoccaFeatureBuilder(file_list, feature_type, final_folder, n_epochs, workers)
 ```

 All three of these processes produce binary files containing the feature data of all channels organized by the subject, session, and epoch. This was done as the COH feature produce a very large number of channels, in excess of the number of epochs that could be produced from the sampled data.
 
## Organize epochs

Adherence to cross validation is carried out within the algorithm evaluation. Data is split into train and test sets that conforms as closely as possible to the initial exhaustive 6 fold cross validation. This is achieved by performing 6 fold cross validation without replacement of the test epochs ensuring each test is unique. As the data is stored by epoch, it is possible to load test and training simultaneously. This process is further streamlined by altering the original data formats into a common format and file structure: **features** X **channels** in each file that is organized as **subject**/**subject_session_epochsxx.bin**.

## Experiments

In the paper two experiments were carried out. **Experiment 1** followed the protocol of La Rocca's testing each individual channel to build an optimal channel subset using their match score-fusion. **Experiment 2** used the PSD features to test the impact of classification as a function of epoch duration using the available trial data.

### Experiment 1

### Experiment 2

```
experimentTrio(param_file,workers)
```

## Results

Performance of the algorithms is evaluated in terms of correct recognition rate (CRR) and equal error rate (EER) averaged over each subject-channel's CV steps. This provides a robust performance metric that can be compared against other EEG classification papers (not all experiments report CCR and/or EER thus having both helps readily compare techniques).

## Major Components

### Mahalanobis Distance

### GMM-UBM

### I-Vector


