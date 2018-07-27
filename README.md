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

To produce the cepstrum (CEP), power spectral density (PSD), and spectral coherence (COH) features the following functions are used to modify the EDF files.

The PSD and COH features were based on the work of La Rocca et al's 2014 paper: _Human Brain Distinctiveness Based on EEG Spectral Coherence Connectivity_ 

 ```
 laRoccaFeatureBuilder(file_list, feature_type, final_folder, n_epochs, workers)
 ```
 
 Given the a **file_list** of EDF files and a desired **feature_type** (0 for PSD, 1 for COH) it will partition the data of each channel into **n_epochs**. They will be stored in **final_folder** (full path name) and be carried out by the number of specified **workers**. While written for the PhysioNet EEGMMIDB, the functionality should extend to any set of EDF files.
 
 To produce the CEP features, an external toolbox is required: [VOICEBOX](http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html). This kit, written by Mike Brookes, provides tools to convert the EDF signals to cepstrum htk files and then read those .HTK files. This is original a tool for speech processing, but why re-invent what already  exists.
 
 To produce the HTK features the same function is used, but **feature_type** is set to 2 for CEP. This approach requires the incoming list to be a list of .HTK files from the already processed .EDF files. Producing the HTK features, 1 second windows with 1/10 second frames using a TCP montage resulting in 22 channels, is outlined in a [paper](https://arxiv.org/ftp/arxiv/papers/1801/1801.02477.pdf) by Harati et al.
 
 ### Produce CEPSTRUM features
 
 The functionality to produce COH and PSD features is built from processing .EDF files. To produce the CEP features, a bit more work is required as the tools to generate them cannot be released at this time. However, to recreate them access to VOICEBOX's melcepst.m, writehtk.m, and readhtk.m functions are all that is required. Using these tools and the associated paper cited in the article allows for the .EDF files to be converted into the features used in experiment. Sample code of how this process is provided below.
```
[edf_header, edf_data] = edfread(file_name);
sample_rate = edf_header.samples(1);
mel_options = 'MtaEfdD';
mel_coef = 8;
sample_end = sample_start + window_size - 1;
data_mel = edf_data(channel,sample_start:sample_end);
htk_features = melcepst(data_mel, sample_rate, mel_options, mel_coef)
epoch_index = [5:10:size(htk_features,1)];
htk_epochs = htk_features(epoch_index,:);
```

 All three of these processes produce binary files containing the feature data of all channels organized by the subject, session, and epoch. This was done as the COH feature produce a very large number of channels, in excess of the number of epochs that could be produced from the sampled data.
 
## Organize epochs

Adherence to cross validation is carried out within the algorithm evaluation. Data is split into train and test sets that conforms as closely as possible to the initial exhaustive 6 fold cross validation. This is achieved by performing 6 fold cross validation without replacement of the test epochs ensuring each test is unique. As the data is stored by epoch, it is possible to load test and training simultaneously. This process is further streamlined by altering the original data formats into a common format and file structure: **features** X **channels** in each file that is organized as **subject**/**subject_session_epochsxx.bin**.

## Experiments

In the paper two experiments were carried out. **Experiment 1** followed the protocol of La Rocca's testing each individual channel to build an optimal channel subset using their match score-fusion. **Experiment 2** used the PSD features to test the impact of classification as a function of epoch duration using the available trial data.

### Experiment 1

```
[CRR, mCRR, EER, mEER, mFPR, mFNR, f_CRR, f_mCRR, CRR_bench] = laRoccaControl(varargin)

where varargin controls Mahal or GMM-UBM/I-Vector analysis

varargin @ Mahala -> file_list, fusion_flag, save_folder, workers
varargin @ GMM-UBM or I-Vector -> file_list, fusion_flag, save_folder, workers, mixtures, iterations, ds_factor, eval_flag
```
This experiment does not run from a parameter file so all arguments must be passed with the function. **file_list** is the full path to a list of the features to be evaluated. **fusion_flag** controls if the match score-fusion function is carried out, where 0 is off and anything else is on. **save_folder** provides the location of where the results will be stored. **workers** tells Matlab how large to make the parallel computing workspace. **mixtures** [2 4 8 ... 512] provides the size of Universal Background Models to build and evaluate. **iteations** controls the number of iterations for producing the UBMs and TVMs. **ds_factor** controls re-sampling when building the UBMs, this is generally set to 1. **eval_flag** controls GMM-UBM (0) or I-Vector evaluation (1).

Each experiment should be saved in a separate file location. The tools used to process the results assume a single experiment per folder in order to produce the plots seen in the paper.
### Experiment 2

```
experimentTrio(param_file,workers)
```

An example paramFile.dat is provided and contains the following:
```
# Parameter File for experiments
#
# subjects
109
#
#
# files, list as many as necessary. assumes they are sessions!
2
C:\Users\NIL\Documents\_NeuroNixCopy\physioNet\_laRoccaFeats\_10second\CEP\features_R01.list
C:\Users\NIL\Documents\_NeuroNixCopy\physioNet\_laRoccaFeats\_10second\CEP\features_R02.list
#
#
# number of iterations to perform experiments
6
#
#
# elements to withhold for testing data
1
#
#
# number of mixtures for UBM
2 4
#
#
# number of iterations to build UBM
20
#
#
# ds_factor, scaling for UBM generation
1
#
#
# save_location, where to write experiment using full pathname
C:\Users\NIL\Documents\_NeuroNixCopy\PhysioNet\paramTest
#
#
#
EOF

```
Each experiment should be saved in a separate file location. The tools used to process the results assume a single experiment per folder in order to produce the plots seen in the paper.

## Results

Performance of the algorithms is evaluated in terms of correct recognition rate (CRR) and equal error rate (EER) averaged over each subject-channel's CV steps. This provides a robust performance metric that can be compared against other EEG classification papers (not all experiments report CCR and/or EER thus having both helps readily compare techniques).

Two functions are used to produce results.
```
crrPlot(condition)
epochVcrr(epoch_10, epoch_5, epoch_2, epoch_1)
```

To use crrPlot() the variable **folders** in the switch statement must be updated to match the proper evaluations. Once mapped, the function operations according to specific figure iterations tried for the paper.

epochVcrr() is more function in that it takes in four experiment folders. It is setup to process epoch duration data from the same series of feature/dataset experiments. Thus each folder should be produced using the same features and base dataset with only a difference in the epoch duration. 

## Major Components

Within these experiments are the major tools to build and evaluate the Mahalanobis Distance classifier from La Rocca, produce Universal Background Models, evaluate GMM-UBMs, and evaluate I-Vectors. The functionality of these tools is limited to the scope of the paper/experiments. However, updated class specific tools (a Universal Background Model class and Total Variability Matrix class) are in progress. These should streamline the operation of UBMs/TVMs/GMM/I-Vector development for the Matlab platform.

### Mahalanobis Distance

```
[CRR, EER, scores, FPR, FNR] = mahalResults(channel_count, epochs, full_index, subject_count, feature_count)
[SUB_SES, SUB, SES, scores] = mahalCcrSlim(training_data, testing_data, subjects, sessions)
```

### Universal Background Models

```
ubms = allUBMs(data_list, nmix, final_niter, ds_factor, folder, verbose)
```

### GMM-UBM

```
[CRR, EER, scores, FPR, FNR] = ubmGmmEvalRocca(training_data, testing_data, ubms, channel)
[SUB_SES, SUB, SES, scores] = ubmGmmEvalRocca2(ubms, training_data, testing_data, subjects, sessions)
```

### I-Vector

```
[CRR, EER, scores, FPR, FNR] = iVectorEvalRocca(training_data, testing_data, ubms,channel)
[SUB_SES, SUB, SES, scores] = iVectorEvalRocca2(ubms, training_data, testing_data, subjects, sessions)
```