# Matlab based GMMUBM and Identity-Vector Classification for EEGs:

The files in this repository are templated from the [MSR Identitiy Toolbox] (http://research.microsoft.com/en-us/downloads/2476c44a-1f63-4fe0-b805-8c2de395bb2c/). After experimenting with the MSR Toolbox I found it lack and began to create something to better suit my needs in terms of ease of use and compliance with newer versions of Matlab.
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

## Evaluate epochs

It is possible to evaluate all features of a given epoch duration against all algorithms.

```
oneLaRocca(exp_folder, list_folder, trial, block_size, covar_flag, workers)
```

The results of all these experiments will be saved into folders created beneath the **exp_folder** directory. The *list_folder* is search for CEP, COH, and PSD folders from which to draw the appropriate feature sets. It is suggested that the feature folder hierarchy is /XX_sec/TYPE which ensures all features are of the same duration.

### Mahalanobis Distance

### GMM-UBM

### I-Vector

## Produce results

Performance of the algorithms is evaluated in terms of correct classification rate (CCR) and equal error rate (EER) averaged over each subject-channel's CV steps. This provides a robust performance metric that can be compared against other EEG classification papers (not all experiments report CCR and/or EER thus having both helps readily compare techniques).

----------
----------

# Legacy README

## Intended Use:

---------------

The data needs to be stored in the *.edf* format with known channel listing. A file list should exist that will point the software to the neccesary files to build all
of the Universal Background Models. These models will all be saved in a user specified folder. The edf data is turned into MelCepst coefficients which is controlled by
the user and stored in a user specified folder. It is encouraged to make separate folders for each run of this system to ensure things are kepy organized and no data is
overwritten. The user will need to know how many cores their computational machine has access to in order to properly run the parallel matlab processes.

## Setup:

---------

These experiment were carried out on a cluster containing 128 cores with 2TB of ram. It is not advisible to run these on a local machine as many facets of the program are
organized to run in parallel. All of the files I contributed were written in Matlab 2016a. However, support from the [Voicebox Toolbox](http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html)
may be helpful in providing auxiliary functions for develop.

## Execution:

------------

The work piggybacked off of concurent research in spike and wave detection for an autoEEG reader. This provided premade [HTK files](http://www.ee.columbia.edu/ln/LabROSA/doc/HTKBook21/node58.html)
previously used as features for that research. If the data type is raw .EDF, then they will need to be read and coverted to Mel Coefficients for use with this software.
It is best to use [edfread](http://www.mathworks.com/matlabcentral/fileexchange/31900-edfread) to process .edf files and then the *melcepst* from Voicebox to convert them
into the desired mel coef format.

### Build Data Matrix

On NeuroNix the software is configured to run off of htk.list files. These were generated by Dr. Picone's group for the AutoEEG development. The script file should look
something like the following:
```
#!/bin/bash
matlab -nodesktop -nosplash -r "buildSubjectSets('/data/isip/exp/tuh_eeg/exp_0166/lists/train_htk.list','_exp/exp_0008/00_raw',32);exit;"
```
This script file should be saved as something, p00_raw.sh, and run with appropriate *qsub* commands. **buildSubjectSets** turns all the raw data into a very large coefficent matrix
of a **.mat** format for use in the ensuing functions. It should be saved to its own directory, see as the second argument passed to the function. The function should be able to make
the directory, but it is best if you create it first. The final argument passed in specifies the number of workers for Matlab's parallelpool to request from the node. This number must
be included in the qsub command or else the problem will either crash or take _forever_ to run. This can also be run directly from Matlab on a local machine, but the number of workers
will probably be closer to 2 or 4. Please be aware that the file you save will be very large as well, and is going to scale depending on the parameters of the files from the list passed.
Please be sure to read up on **buildSubjectSets** below for more information.

### Generate Universal Background Models

Once the list of files is converted into the data matrix, the *universal background models* must be generated from the data. These seed any and all possible discrimintation of the data.
Acting as the basis of both *Joint Factor Analysis* and *Identity Vectors*, allows for more variation in processing the data. A typical command should look like the follow:
```
#!/bin/bash
matlab -nodesktop -nosplash -r "a=FindMatchingFiles('_exp/exp_0008/00_raw','mel');load(a{1});htkUbmGen(mel_data,[2 4 8 16 32 64 128 256],'_exp/exp_0008/01_ubm',32,10,1);exit;"
```
This could be streamlined if you wish to hardcode in the filename to load. The second argument passed in are the various Gaussian mixtures to be modeled, powers of 2 only. Then the next two
are paramteres for the search in terms of iteration number, 10, of the model generating and a subsampling factor for the expectation maximization tool that creates the GMM. For each parameter
in the original .htk files, an n-dimensional Gaussian is generated to cover the entire dimensional space. For 2 mixtures, only two Gaussians are built that must cover/classify all the presented
data. In this manner, the data could be overfit (and under fit) easily if attention is not paid to the original feature count.

At this point, the UBM results will indicate error rates across the data set in terms of Equal Error Rate. This value shows the corresponding flase positive and false negative rates which is
somewhat helpful in understanding the GMMs accuracy. Ideally, an error rate for specific models should be generated against each data set. **I haven't been able to write this yet.**

### Build Identity Vectors

The final action is to generate Identity Vectors from the data. The command to achieve this looks as follows:
```
#!/bin/bash
matlab -nodesktop -nosplash -r "feedIvector('_exp/exp_0008/01_ubm,'_exp/exp_0008/00_raw',32,'_exp/exp_0008/02_iVec');exit;"
```
As the present final command, this generates I-Vectors and evaluates them against the entire training set. Since the goal is only to match, there is no need to split test and training as the
program should never be asked to find specific segments. We're trying to match large segements of recordings, but perhaps smaller sections could work as well? Again, this needs the directories
of all the data along with a worker limit for Matlab's parallel processing tools. The results come out as EERs and a .csv table of false positives and true positives. A very basic ranking system
has been implemented that attempts to determine how _alike_ all the subjects are to each other. As this entire process treats channels as new idenpendent entries there is work to be done on
understanding how channel dependency impacts the results.

At this point you've run all the major functions I have written and adapted for my research. Ideally this should all work, if not please let me know so I can fix whatever isn't working.

## File Listing:

----------------

- [buildSjubectSets](#buildSubjectSets): Works
- [calculateUBMerrors](#calculateUBMerrors): Works, obsolete
- [channelMatch](#channelMatch): Broken, repair
- [compute_eer_2](#compute_eer_2): Works, not written by me
- [cosineDistance](#cosineDistance): Works
- [directoryCheck](#directoryCheck): Works
- [edfFolderGen](#edfFolderGen): Works, obsolete
- [edfFolderGenWinFra](#edfFolderGenWinFra): Works
- [edfread](#edfread): Works, not written by me
- [enhanced_rdir](#enhanced_rdir): Works, not written by me
- [errorPlot](#errorPlot): Works, obsolete
- [evaluateIvector](#evaluateIvector): Works
- [feedIvector](#feedIvector): Works
- [fileListConvert](#fileListConvert): Works
- [fileNameMatch](#fileNameMatch): Works
- [findDirectoryMatch](#findDirectoryMatch): Works
- [findIvectorMatch](#findIvectorMatch): Works
- [findMatchingFiles](#findMatchingFiles): Works
- [folderCount](#folderCount): Works, obsolete
- [folderFinder](#folderFinder): Works
- [folderNameCount](#folderNameCount): Works
- [folderTestSplit](#folderTestSplit): Works, obsolete
- [fullEDFbuild](#fullEDFbuild): Works, obsolete
- [fullSubjectGroupTrain](#fullSubjectGroupTrain): Works, obsolete
- [fullSubjectTrain](#fullSubjectTrain): Works, obsolete
- [fullUbmBuild](#fullUbmBuild): Works, obsolete
- [genUbmData](#genUbmData): Works, obsolete
- [getAllFiles](#getAllFiles): Works, not written by me
- [gmm_demo](#gmm_demo): Works, obsolete
- [htkToUbm](#htkToUbm): Works
- [htkToUbmOnly](#htkToUbmOnly): Works
- [htkUbmGen](#htkUbmGen): Works
- [htkUbmGenOnly](#htkUbmGenOnly): Works
- [iVectorWCompare](#iVectorWCompare): Works
- [makeUBM](#makeUBM): Works
- [makeUBMnoPlots](#makeUBMnoPlots): Works
- [makeUBMOnly](#makeUBMOnly): Works
- [melTestTrainSplit](#melTestTrainSplit): Works
- [melTestTrainSplit_2](#melTestTrainSplit_2): Works
- [modelAccuracy](#modelAccuracy): Works
- [poopulateUBMFolder](#poopulateUBMFolder): Works
- [processDirectory](#processDirectory): Works
- [processEDF](#processEDF): Works
- [processEDFsingle](#processEDFsingle): Works
- [processEdfWindowFrame](#processEdfWindowFrame): Works
- [processsIvector](#processsIvector): Works
- [processSPH](#processSPH): Works
- [processTestTrainUbm](#processTestTrainUbm): Works
- [processUBM](#processUBM): Works
- [rankedSubjects](#rankedSubjects): Works
- [rdir](#rdir): Works, not written by me
- [scoreIvector](#scoreIvector)
- [singleChannelEDF](#singleChannelEDF)
- [splitFeatureFolders](#splitFeatureFolders)
- [splitFeatureWinFram](#splitFeatureWinFram)
- [splitTrainTest](#splitTrainTest)
- [spiltUBMfolders](#spiltUBMfolders)
- [testErrorPlot](#testErrorPlot)
- [testThisThing](#testThisThing)
- [testUBM](#testUBM)

### <a name="buildSubjectSets">buildSubjectSets</a>
```
mel_data = buildSubjectSets(file_list, save_folder, workers)
```
- **mel_data**: resultant matrix of .htk files
- **file_list**: the path to the file of interest, should be an universal path name not a local path name
- **save_folder**: location where *mel_data* will be saved
- **workers**: number of workers allowed during parfor loops

Function created to parse ISIP file lists containing locations of .htk files. Originally the program operated by passing in a folder name and
using all .edf files found within that folder. While this worked for local testing, such a system was not feasible for NeuroNix. Instead I
switched over to using the same 'system' employed by the ISIP group which relies on a file containing the location of all the pertinent data.
This function reads the given file and generates the necessary data matrix for later steps of the process.

### <a name="calculateUBMerrors">calculateUBMerrors</a>
```
error_results = calculateUBMerrors(folder_name)
```
- **error_results**: returns all the error rates, EERs, in a matrix
- **folder_name**: pass in a folder containing the test and training data to be evaluated along with the ubm data, searches recursively in the
folder directory

This is essentially legacy code from initial testing on the physionet motor database. Everything was segemented into test/training so I had to
find some way to evaluate it locally. It was all stored in folders and I made separate experiment folders for each test which held the test and
training data for the subjects.

### <a name="channelMatch">channelMatch</a>
```
diff_weighted = channelMatch(gmm,channel_data)
```
- **diff_weighted**: distance between the gmm and the channel data
- **gmm**: structure containing a generated UBM
- **channel_data**: data for a specific channel/subject to be compared against the UBM

The concept was to understand how the UBMs were capturing the channel data in an effort to build a sample by sample heatmap of similiarity.
Ideally finding matches is only useful if those matches can teach you something about what is being compared. It seemed reasonable to build a
sample by sample heat map that may allow insight when overlapped with AutoEEG. Unfortunately this doesn't work because the structure of the data
and the UBM don't lend to a straightfoward calculation. I thought they would, but I must have missed something; I need to look at this again.

### <a name="compute_eer_2">compute_eer_2</a>
```
[eer, dcf08, dcf10] = compute_eer_2(scores, labels, showfig)
```
- **eer**: percent equal error rate
- **dcf08**: minimum detection cost function with SRE'08 parameters
- **dcf10**: minimum DCF with SRE'10 parameters
- **scores**: likelihood scores generated between target and non-target trials
- **labels**: index of what matches to what for evaluation purposes, 0 & 10
- **showfig**: flag to display EER plot

This is a modified version of *computer_eer* in only the plot function. The original would only plot errors below 40% whereas this plots errors
up to 0.95 percent. It was necessary for initial testing as the errors were quite large when I was learning the process. The original function,
*computer_eer*, is part of the Microsoft I-Vector package.

### <a name="cosineDistance">cosineDistance</a>
```
result = cosineDistance(i_vector)
```
- **result**: distances between the I-Vectors
- **i_vector**: array of I-Vectors to compared

Based upon the literature, the easiest way to compare I-Vectors is the cosine of the angle between their vectors. This amounts to the dot products
and normal products producing a value fed into arccos.If the resulant angle is large they don't point in the same directions, but smaller angles
suggest the vectors are similar. Ideally this should show subjects recording in the same environment being alike, but could also provide insight
into what types of featres are most dominant in making matches OR if any features are linked. This could be an interesting area of research, but
would rely heavily on Joint Factor Analsys to separate the data into eigenvoice, eigenchannel and residual data segments.

### <a name="directoryCheck">directoryCheck</a>
```
directoryCheck(directory_path)
```
- **directory_path**: targeted save path for a file

Checks to see if the directory is valid. If it is not, it creates a directory and alerts the user via print statements.

### <a name="edfFolderGen">edfFolderGen</a>
```
edfFolderGen(file_names, target_folder, channels, window_size, mel_coefs)
```
- **file_names**: list of all the .edf file names to be processEDF
- **target_folder**: location of where data will be written
- **channels**: number of channels to process in each .edf file
- **window_size**: variable for window size used in the mel coefficient generator
- **mel_coefs**: variable to specify the number of mel coefficients to generate

Original core function that processed folders of .edf files. Used a static frame size and only allowed the user to control the *window_size* and
*mel_coef* count. Would only generate mel coefficient data and save it in the desired folder. Minimal control over parameters but used to quickly
setup a test system on local data.

### <a name="edfFolderGenWinFra">edfFolderGenWinFra</a>
```
edfFolderGenwinFra(file_names, target_folder, mel_coefs, window_size, frame_size, workers)
```
- **file_names**: list of all the .edf file names to be processEDF
- **target_folder**: location of where data will be written
- **mel_coefs**: variable to specify the number of mel coefficients to generate
- **window_size**: variable for window size used in the mel coefficient generator
- **frame_size**: variable to specify the size of the processing frame
- **workers**: variable to enable parfor loops within Matlab

Updated version of *edfFolderGen* that now enables proper frame/window handling of incoming data. Still only able to read .edf files and creates
similar output of mel coefficient matrix as older function.

### <a name="edfread">edfread</a>
```
[hdr, record] = edfread(fname, varargin)
```
- **hdr**: header file from edf
- **record**: raw data from edf
- **fname**: file name of edf
- **varargin**: variable parameters to pass to function

Enables the European Data Format files to be parased by Matlab into headers and raw data.

### <a name="enhanced_rdir">enhanced_rdir</a>
```
d = rdir(path)
```
- **path**: either blank or specified folder path

A third party enhanced version of Matlab's built in dir() function.

### <a name="errorPlot">errorPlot</a>
```
errorPlot(folder_name)
```
- **folder_name**: location of output files for debugging

Plots error rates of various stages of the process after completion. Should show error from PLDA, UBM and I-vectors assuming said data is generated
by higher level functions.

### <a name="evaluateIvector">evaluateIvector</a>
```
[eer, model_IVs, final_test_IVs] = evaluateIvector(ubm_data, features, test_split, train_split, workers)
```
- **eer**: equal error rate result
- **model_VIs**: the chosen model I-Vectors for each of the subjects
- **final_test_IVs**: the I-Vectors used to test against the model I-Vectors
- **ubm_data**: a single structure of the UBM pertaining to one mixture variant
- **features**: the raw feature data matrix
- **test_split**: the percentage of the data that should be turned into testing data
- **train_split**: the percentage of the data that should be turned into training data
- **workers**: the number of workes allocated for parfor loops 

Generated I-Vectors from the raw data and tested them against chunks of testing data. A very poor initial function as varying test/training sizes is
not wise.

### <a name="feedIvector">feedIvector</a>
```
feedIvector(ubm_folder, mel_folder, workers, save_folder)
```
- **ubm_folder**: path to the folder containing all the UBM data
- **mel_folder**: path to the folder containing the raw mel coef data
- **workers**: number of workers allowed for in parfor loops
- **save_folder**: path to the folder where results will be written

This function assumes htk files are being targeted and is meant to operate on NeuroNix. It needs a folder structure as shown in the *Execution* and is
ideally suited for cluster operation. It writes a number of files to the result directory through the *scoreIvector* function and attempts to handle
naming conventions with the internal function *mixtureSplit*, `result=mixtureSplit(file_list)`, that chomps through the file names to ensure everything
is saved in folders sorted by the UBM mixture count.

### <a name="fileListConvert">fileListConvert</a>
```
[result, count] = fileListConvert(file_name)
```
- **result**: the returned value from textscan() on the given file
- **count**: returns the number of file listings
- **file_name**: path to the file containing all the locations for the htk list files

Reads path names from a list file and gathers them in a cell array for processing via the existing software tools

### <a name="fileNameMatch">fileNameMatch</a>
```
new_path_ways = fileNameMatch(file_paths, queries)
```
- **new_path_ways**: numerical portion of the file paths passed into the function
- **file_paths**: array of file paths to search
- **queries**: array of size groups by trials for searching the PhysioNet data

Written to support processing data from the PhysioNet EDFs. Designed around ordering the subject/trial folder hierarchy.

### <a name="findDirectoryMatch">findDirectoryMatch</a>
```
result = findDirectoryMatch(query, varargin)
```
- **result**: a cell of folders matching the search queries
- **query**: phrasing to match in search
- **varargin**: if nothing, uses the present directory otherwise assumes you've handed it a directory to search

A helper function that returns only file names matching a given criteria.

### <a name="findIvectorMatch">findIvectorMatch</a>
```
[pruned_result, t_p_percent] = findIvectorMatch(train, test)
```
- **pruned_result**: Accuracy of I-Vector against the test set
- **t_p_percent**: the true positive percent of each of the given I-Vectors

A tool to evaluate I-Vector performance of the 'chosen' test markers against all of the generated training vectors.

### <a name="findMatchingFiles">findMatchingFiles</a>
```
data_files = findMatchingFiles(folder_name, varargin)
```
- **data_files**: file list of all files matching search query
- **folder_name**: folder to search through
- **varargin**: arguments to be checked for in the search

Each additional input is searched for through the results of the recursive file search. The varargin are used to prune the results starting with the
1st entry.

### <a name="folderCount">folderCount</a>
```
folder_count = folderCount(input_string)
```
- **folder_count**: number of folders matching input_string in the present directory
- **input_string**: string to match to folder names

Counts how many folders match the specified input_string in the present directory.

### <a name="folderFinder">folderFinder</a>
```
data_folders = folderFinder(folder_name, key_word, varargin)
```
- **data_folders**: a list of all folders matching the search criteria
- **folder_name**: the folder where the search should starting
- **varargin**: just a flag to make the search recursive in the specified folder or not

Searches folders and returns only those matching the keyword. Able to search recursively or just one level deep. 

### <a name="folderNameCount">folderNameCount</a>
```
result = folderNameCount(input_string, folder_name)
```
- **result**: number of folders matching the input string
- **input_string**: query to match to folder name
- **folder_name**: directory to start folder search within

A variation on the previously mentioned *folderCount*, but now the user can specify the folder to search.

### <a name="folderTestSplit">folderTestSplit</a>
```
folderTestSplit(folder_name, split, file_name)
```
- **folder_name**: location of the *data_* file
- **split**: percentage split to be made testing data
- **file_name**: file name of original data

Fragment of early code used to test on the early local data sets.

### <a name="fullEDFbuild">fullEDFbuild</a>
```
fullEDFbuild(file_name, channels, segements, mel_win, mel_coef, split_percent, workers)
```
- **file_name**: path to file
- **channels**: number of channels present in the file
- **segments**: number of segements to break the file into (this was a bad idea)
- **mel_win**: size of the mel coefficient window in milliseconds
- **mel_coef**: number of coefficients to make
- **split_percent**: how much testing dating to generate
- **workers**: limit for parallel functions in matlab

Original function to create the coefficient matrix.

### <a name="fullSubjectGroupTrain">fullSubjectGroupTrain</a>
```
fullSubjectGroupTrain(root_directory, split, grouping)
```
- **root_directory**: directory containing all the processed files
- **split**: testing split percentage
- **grouping**: key word to organize the data into know sections, based upon subject?

Orignal function to train large file sets after the coef matrix was built.

### <a name="fullSubjectTrain">fullSubjectTrain</a>
```
fullSubjectTrain(root_directory, split)
```
- **root_directory**: location of the processed data
- **split**: percentage split of the testing data

Combines all of the subject data, which was once done trial by trial, to develop a set to train across all of the subect's trials.

### <a name="fullUbmBuild">fullUbmBuild</a>
```
fullUbmBuild(split_folder, mixtures, split, workers)
```
- **split_folder**: folder containing the separated data
- **mixtures**: array of n-GMM desired mixtures
- **split**: how the split value is for filename generation
- **workers**: how many workers can be assigned for parallel processes

Would take a folder with test and train data then produce *n* UBMs in the specified results folder(internally set name).
 
### <a name="genUbmData">genUbmData</a>
```
genUbmData(feature_folder, test_split, train_split, mixtures, channels, workers)
```
- **feature_folder**: specify location of feature data
- **test_split**: percentage of data to be turned into testing
- **train_split**: percentage of data to be turned into training
- **mixtures**: specify the desired Guassian mixtures
- **channels**: indicate how many channels should be used
- **workers**: allocate parallel process worker number

Function to handle calling and saving results from generate UBMs. Mainly used with initial local testing of EDF data.

### <a name="getAllFiles">getAllFiles</a>
```
fileList = getAllFiles(dirName)
```
- **fileList**: a list of all files within the directory
- **dirName**: directory to scrape for files

Returns all files, recursively, within a given directory.

### <a name="gmm_demo">gmm_demo</a>
```
gmm_demo
```

A modified version of the GMM-UBM demo from the Microsoft Identity Toolbox. Used mostly to understand the links between the individual functions.

### <a name="htkToUbm">htkToUbm</a>
```
htkToUbm(file_list, workers, mixtures, save_folder, iterations, ds_factor)
```
- **file_list**: the htk.list file containg all the locations of the processed data
- **workers**: number of parallel operations to allow
- **mixtures**: array contain the size of the desired GMM-UBMs
- **save_folder**: location to save results of UBMs
- **iterations**: number of times to optimize the UBM
- **ds_factor**: scaling factor to allow for sliding between samples

A helper function designed to break up commands in the shell script on NeuroNix. Would build the subject data matrix and then generate associated UBMs.


### <a name="htkToUbmOnly">htkToUbmOnly</a>
```
htkToUbmOnly(folder_name, name, train_data, mixture, iterations, ds_factor, workers)
```
- **folder_name**: the htk.list file containg all the locations of the processed data
- **name**: for saving the approriate file
- **train_data**: the targeted data to train the models on
- **mixture**: the desired GMM-UBMs
- **iterations**: number of times to optimize the UBM
- **ds_factor**: scaling factor to allow for sliding between samples
- **workers**: number of parallel operations to allow

A helper function that only calls *gmm_em* to build the a single UBM for the passed in mixture.

### <a name="htkUbmGen">htkUbmGen</a>
```
htkUbmGen(coefs, mixtures, save_folder, workers, iterations, ds_factor)
```
- **coefs**: coef variable matrix
- **mixtures**: array containing the UBMs to generate
- **save_folder**: location to save the UBMs
- **workers**: number of parallel options to enable
- **iterations**: number of times to process the final UBM
- **ds_factor**: scaling factor for how to slide sample during processing

The main function that sets up and saves the results from processTestTrainUbm. Also records error results.

### <a name="htkUbmGenOnly">htkUbmGenOnly</a>
```
htkUbmGen(coefs, mixtures, save_folder, workers, iterations, ds_factor)
```
- **coefs**: coef variable matrix
- **mixtures**: array containing the UBMs to generate
- **save_folder**: location to save the UBMs
- **workers**: number of parallel options to enable
- **iterations**: number of times to process the final UBM
- **ds_factor**: scaling factor for how to slide sample during processing

The main function that sets up and saves the results from processTestTrainUbm. Does not record error results, only saves UBM.

### <a name="iVectorWCompare">iVectorWCompare</a>
```
[result, ttest_8, ttest_16] = iVectorWCompare()
```
- **result**: distance between I-Vectors
- **ttest_8**: T test result of I-Vectors
- **ttest_16**: T test result of I-Vectors

A minor tool to compute the difference between I-Vectors in a folder. Hardcoded to find directory matching 'single'. Then plots results of all I-Vectors
found.	

### <a name="makeUBM">makeUBM</a>
```
[ubm, eer] = makeUBM(folder_name, train_file, test_file, mixtures, iterations, ds_factor, workers)
```
- **folder_name**: location to save UBMs
- **train_file**: file containing training data
- **test_file**: file containing testing data
- **mixtures**: array containing the UBMs to generate
- **iterations**: number of times to process the final UBM
- **ds_factor**: scaling factor for how to slide sample during processing
- **workers**: number of parallel options to enable

Original function to produce UBMs once training and testing data had been split from raw coefficient data. Produces a lot of plots which triggered a minor
memory leak when run on NeuroNix. Produces error results from testing data.

### <a name="makeUBMnoPlots">makeUBMnoPlots</a>
```
[ubm, eer] = makeUBMnoPlots(folder_name, train_file, test_file, mixtures, iterations, ds_factor, workers)
```
- **folder_name**: location to save UBMs
- **train_file**: file containing training data
- **test_file**: file containing testing data
- **mixtures**: array containing the UBMs to generate
- **iterations**: number of times to process the final UBM
- **ds_factor**: scaling factor for how to slide sample during processing
- **workers**: number of parallel options to enable

Original function to produce UBMs once training and testing data had been split from raw coefficient data. Produces error results from testing data.

### <a name="makeUBMOnly">makeUBMOnly</a>
```
ubm = makeUBMonly(folder_name, train_file,  mixtures, iterations, ds_factor, workers)
```
- **folder_name**: location to save UBMs
- **train_file**: file containing training data
- **mixtures**: array containing the UBMs to generate
- **iterations**: number of times to process the final UBM
- **ds_factor**: scaling factor for how to slide sample during processing
- **workers**: number of parallel options to enable

Stand alone function to produce results for one UBM given only training data. Does not evaulate results.

### <a name="melTestTrainSplit">melTestTrainSplit</a>
```
[test_set,train_set] = melTestTrainSplit( coefficients, split_percent, data_folder, file_name )
```
- **test_set**: variable holding all the test data
- **train_set**: variable holging all the train data
- **coefficients**: variable holding all of raw mel coef data
- **split_percent**: amount of data to use for test
- **data_folder**: location to save splits
- **file_name**: original file name of data being analyzed

Splits data up according to the segements by removing whole chunks. The result is a non-random sampling of the data.

### <a name="melTestTrainSplit_2">melTestTrainSplit_2</a>
```
[test_set,train_set] = melTestTrainSplit2( coefficients, split_percent, data_folder, file_name )
```
- **test_set**: variable holding all the test data
- **train_set**: variable holging all the train data
- **coefficients**: variable holding all of raw mel coef data
- **split_percent**: amount of data to use for test
- **data_folder**: location to save splits
- **file_name**: original file name of data being analyzed

Rebuilds data out of segments and randomly pulls individual samples for testing.

### <a name="modelAccuracy">modelAccuracy</a>
```
result = modelAccuracy(model_results)
```
- **result**: an n-subject by n-subject array showing the hit/miss of each model
- **model_results**: array holding the scores of each subject/sample combination

Produces a more granular statistic than the EER associatd with both the UBM model and the I-Vector model.

### <a name="poopulateUBMFolder">poopulateUBMFolder</a>
```
[ubm_results, error_ubm]=populateUBMfolder(folder_name, workers, full, varargin)
```
- **ubm_results**: the UBMs for each mixture
- **error_ubm**: the errors associated with each UBM
- **folder_name**: location to save results
- **workers**: number of parallel works to ask for in parfor loops
- **varargin**: an optional input to specify the n-mixtures, defaults to [2,4,8,16,32,64,128,256]

A function to allow folders full of test/train data to be quickly turned over into UBM data. Initially used locally, but then also deployed on NeuroNix.

### <a name="processDirectory">processDirectory</a>
```
[error_PLDA, true_positive, iv_map] = processDirectory(folder_name, varargin)
```
- **error_PLDA**: error values from the I-Vector analysis
- **true_positive**: true positives associated with each I-Vector model
- **iv_map**: I-Vector mapping used for analysis
- **folder_name**: folder to use to find UBMs and begin the analysis
- **varargin**: allows the number of workers to be specified, otherwise defaults to 4

Processing for locally building I-Vectors and their results.

### <a name="processEDF">processEDF</a>
```
set_c = processEDF(file_name, channels, segments, mel_window, mel_coef)
```
- **set_c**: the matrix of raw coefs
- **file_name**: the edf file to read
- **channels**: number of channels valid in the file
- **segments**: how many sections to break the file into
- **mel_window**: how large of a sampling window to use, in milliseconds
- **mel_coef**: how many coefficients to use

The original function that would convert EDF into a data matrix for testing.

### <a name="processEDFsingle">processEDFsingle</a>
```
set_c = processEDFSingle(file_name, target_directory, channels, segments, mel_window, mel_coef)
```
- **set_c**: the matrix of raw coefs
- **file_name**: the edf file to read
- **target_directory**: location to save the resultant *set_c* matrix 
- **channels**: number of channels valid in the file
- **segments**: how many sections to break the file into
- **mel_window**: how large of a sampling window to use, in milliseconds
- **mel_coef**: how many coefficients to use

### <a name="processEdfWindowFrame">processEdfWindowFrame</a>
```
result = processEdfWindowFrame(input_data, save_location, rate_samples, size_window, size_frame, count_mel_coef, workers)
```
- **result**: matrix of coef data
- **input_data**: raw data from edf file
- **save_location**: wehre to save the matrix coef data
- **rate_samples**: specify the sampling rate
- **size_windows**: specifiy the window size for sampling
- **size_frame**: specifiy the frame size for sampling
- **count_mel_coef**: number of mel coef to use
- **workers**: number of workers to use in parfor loops

The first function to properly process the raw data by using both windows and frames.

### <a name="processsIvector">processsIvector</a>
```
processIVector(mel_index, ubm_index, tvDim, iterations, workers)
```
- **mel_index**: indicate how many coef were used
- **ubm_index**: indicate with UBM to use via its tag
- **tvDim**: specifiy the dimension of the tv matrix trained
- **iterations**: set the number of iterations for training
- **workers**: enable some number of parfor workers

This is a modified version of the I-Vector generator found in the Microsoft suite of software tools. I mostly used it for testing
and setting up other functions. Uses the *mel_index* and *ubm_index* to build a folder to save the results.

### <a name="processSPH">processSPH</a>
```
c_save_file = processSPH(folder,varargin)
```
- **c_save_file**: output data matrix of features
- **folder**: folder holding sph data
- **varargin**: allows for channels, mel_coef, and mel_window to be specified when generating feature matrix

This allowed for testing on native speech data in the .SPH (sphere) format from the switchboard data set.

### <a name="processTestTrainUbm">processTestTrainUbm</a>
```
eer = processTestTrainUbm(folder_name, name, train_data, test_data,  mixtures, iterations, ds_factor, workers)
```
- **eer*: Equal Error Rate results from trained UBMs
- **folder_name**: name of folder to pull data from
- **name**: name of original data file
- **train_data**: raw data for training in matrix format
- **test_data**: raw data for testing in matrix format
- **mixtures**: vector of GMM-UBM mixtures to generate
- **iterations**: number of times to repeat the UBM algorithm
- **ds_factor**: sliding sample factor for UBM algorithm
- **workers**: number of cpus to allow Matlab to use for parfor loops

This function evaluated a UBM once it was generated and saved the results in a corresponding file.

### <a name="processUBM">processUBM</a>
```
processUBM(mel_index, mixtures, iterations, ds_factor, workers)
```
- **mel_index**: numeric value of what version of the coef data to use
- **mixtures**: vector of desired GMM-UBMs to build
- **iterations**: number of times to repeat the UBM algorithm
- **ds_factor**: sliding sample factor for UBM algorithm
- **workers**: number of cpus to allow Matlab to use for parfor loops

Computes the desired UBMs and evaluates them. The results are saved in a folder generated by the function relating to the mel coefficients used
and the associated GMM-UBM mixture.

### <a name="rankedSubjects">rankedSubjects</a>
```
result = rankedSubjects(data)
```
- **result**: subject-by-subject matrix showing the hit/miss rate of the resultant I-Vector distances
- **data**: a data matrix containing n-subjects by m-samples showing which the order the I-Vector matches

A function to produce a matrix of data that shows weighted averages based upon the distance each I-Vector is from a given sample data set. This should
allow more insight into which I-Vectors are closer to others in an effort to start grouping them for further analysis.

### <a name="rdir">rdir</a>
```
[varargout] = rdir(rootdir,varargin)
```
- **varargout**: generally a file list, but can be two argumetns
- **rootdir**: directory to search
- **varargin**: most often a search time, wildcards accepted

Helper funtion from the online Mathworks Repo that is more powerful than the native *dir* function.

### <a name="scoreIvector">scoreIvector</a>
```
[eer, model_acc, iv_scores, model_IVs, rankedResults ] = scoreIvector(ubm_data, train_data, test_data, workers, save_path)
```
- **eer**: Equal Error Rate results
- **model_acc**: Matrix showing the accuracy of the modeled I-Vectors
- **iv_scores**: Raw I-Vectors scores
- **model_IVs**: The I-Vector models used
- **rankedResults**: ranked results based upon the I-vector scores
- **ubm_data**: variable containing the UBM data
- **train_data**: variable containing the training data
- **test_data**: variable containing the test data (which should be identical to the training data)
- **workers**: number of workers Matlab to call for parfor loops
- **save_path**: folder to save results into

The final function in the process that determines the I-Vector scores and presents all the resultant data.

### <a name="singleChannelEDF">singleChannelEDF</a>
```
[ubm_dat, ubm_err] = singleChannelEDF(folder_name, channel_index, window)
```

- **ubm_dat**: resultant UBMs in struct format
- **ubm_eer**: evaulation of UBM on test data
- **folder_name**: parent directory containing the raw data sets
- **channel_index**: which channels to process
- **window**: size of window in coefficient processing

Designd to enable specific processing of files based upon various channel configurations. Turned out not to be an issue, but brings up interseting
ideas about how this process differs with EEGS than speechs signals.

### <a name="splitFeatureFolders">splitFeatureFolders</a>
```
splitFeaturesFolders(root_directory, split, workers)
```
- **root_directory**: location of data sets to be split
- **split**: percentage of files to turn into test data
- **workers**: number of workers to allocate for parpool.

A quick fix to generating test and train sets after the raw data is processed. Turned out to be useless since data shouldn't be split for speaker recognition
problem.

### <a name="splitFeatureWinFram">splitFeatureWinFram</a>
```
splitFeaturesWinFram(root_directory, split, workers)
```
- **root_directory**: location of data sets to be split
- **split**: percentage of files to turn into test data
- **workers**: number of workers to allocate for parpool.

Splits test and train sets from raw data, but takes into account window and frame sizing.

### <a name="splitTrainTest">splitTrainTest</a>
```
[test_set,training_data] = splitTrainTest( coefficents, train_split, folder, varargin )
```
- **test_set**: cell variable holding test data
- **training_set**: cell variable holding train data
- **coefficients**: raw coefficient data
- **train_split**: how much training data to split, numeric value
- **folder**: save location for cell test/train variables	
- **varargin**: how much training data to make if only one additional argument passed, a fifth argument defaults to half and half

A more complex tool to split test and train sets. Yet again, turns out to be useless once the project gained focus.

### <a name="spiltUBMfolders">spiltUBMfolders</a>
```
splitUBMfolders(file_list,split,mixtures,iterations,ds_factor,workers)
```

**Broken function never completed.**

### <a name="testErrorPlot">testErrorPlot</a>
```
testErrorPlot(test_error,sub)
```
- **test_error**: matrix of erros assocaited with UBMs for a given subject or set of subjects
- **sub**: either 't' or 's' for subjects or trials. ensures the correct labels are generated on the plots

Function written to produce plot for the EMBC 2016 conference paper.

### <a name="testThisThing">testThisThing</a>
```
eer = testThisThing(input_data, ubm)
```

**One time use debug function. Ignore**

### <a name="testUBM">testUBM</a>
```
[eer, fig_handle] = testUBM(train_data,test_data,ubm)
```
- **eer**: Equal Error Rate
- **fig_handle**: figure handle for manipulation
- **train_data**: train data, must match test data
- **test_data**: test data, must match train data
- **ubm**: ubm struct

Function to evaluate UBMs on different subject/channel data. Looking to find differences across subjects/channels, if any.

### <a name="trialDetails>trialDetails</a>
```
[coef,score,latent,tsq,expla,mu] = trialDetails(full_errors)
```

**One time use debugging function. Please ignore.**

### <a name="ttestFullArray">ttestFullArray</a>
```
result = ttestFullArray(data)
```

**One time use debugging function. Please ignore.**

### <a name="writeCellToFile>writeCellToFile</a>
```
writeCellToFile(cell_string,file_name,varargin)
```
- **cell_string**: cell of string arrays to be written
- **file_name**: desired file name
- **varargin**: controls for how to format the file name, five cases	

Allow a cell array of strings to be written to a file.