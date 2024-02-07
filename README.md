# Chronic-Pain-Project

The purpose of this repository is to share pre-processing and data analysis scripts with other researchers. Our study examined predictive coding differences between chronic pain patients and healthy controls. We performed time-domain, time-frequency and source localisation analyses on EEG signals which were collected during an auditory oddball paradigm.

If you have any questions about the scripts, please feel free to contact colum.ose28@gmail.com. I'll do my best to get back to you quickly if you enter something along the lines of 'github script queries' in the subject.

Information pertaining to specific aspects of the processing scripts can be found below:

## Pre-processing
I followed a pretty traditional method of EEG pre-processing. Artifact rejection was done manually in EEGLAB, as automatic methods such as those in Field Trip have provided me with noisy data in the past.

## Time-domain analysis
### Creating sample means
Credit goes out to my friend and colleague Feifan Chen for this!

With conditioning experiments, conditions may be presented with highly uneven ratios, to prime participants to predict
the sensation of a given condition. The downside of this is that comparing average ERPs of conditions with a highly uneven trial count is tricky. The condition with more trials will have less variance in the average ERP.

To correct for this, we created a distribution of sample means, and then used the mean of this distribution as the average ERP for each participant. 

For example...Condition 1 has 200 trials and condition 2 has 50 trials.

Select 50 trials randomly from condition 1, to match the trial count of condition 2. Average these trials to create an average ERP. Repeat this procedure 1,000 times to create a distribution of these permuted means. Then we take the average of the distribution to compare with condition 2.

### Field Trip permutation tests
We used cluster-based permutation paired/independent samples t-tests to detect within and between group differences across all channels and time-points. The benefit of such a method is that it corrects for the family-wise error rate that accumulates across multiple tests. You can read plenty about this from [Field Trip tutorials](https://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/).

#### Plotting
Significant channel-time points were obtained using cluster-based permutation tests. The output of these [scripts](https://github.com/columose/Chronic-Pain-Project/tree/3470c2100d5b2961176fd0bac76008701f737c3f/Plotting%20time%20domain%20results) is a nice topoplot using EEGLAB functions, where significant channels are highlighted in **bold**. The script is nice as it does all the plotting and highlighting dynamically. We also plot a time-series, averaged over the significant channels, where the significant timeframe is shaded in green.

### Source localisation
I performed source localisation of ERPs using [Brainstorm toolbox](https://neuroimage.usc.edu/brainstorm/). The toolbox has a GUI, but also provides you with .mat script after performing operations so you can also code the analyses.

I'd definitely recommend **saving the links** to files in a folder so you can reload them independently of the GUI.

It's a useful toolbox for people, but playing around with the GUI can be time-consuming.

### MNE time-domain analysis
I replicated my Field Trip permutation tests in [MNE](https://mne.tools/stable/index.html) to get a feel for the toolbox. I would try the toolbox again but have a preference for Field Trip as it's more widely used and has superior documentation.


## Time-frequency analysis
## Signal decomposition
I performed time-frequency decomposition of data from a single channel, 'Cz' using an adapted Mike Cohen script.
The adapted version is nice because it deals with sharp edges during trial concatenation using 'reflection'.
I also updated the function that I used to perform a transformation of complex numbers to real numbers, 
a nice tweak that decreases file size and computation time dramatically. 

## Pixel statistics
I performed permutation testing to identify statistically significant time-frequency pixels and then outlined
those pixels on the spectrogram. Another nice idea from [Mike X Cohen!](https://www.youtube.com/watch?v=fAYFtpKwJRQ&list=PLn0OLiymPak1Ch2ce47MqwpIw0x3m6iZ7&index=6)
### Plotting
The [script](https://github.com/columose/Chronic-Pain-Project/blob/e5941e5741be29e28d199885c0165fa33d399a81/Time-frequency%20pixel%20statistics%20and%20plotting/TF_phase_Cz_EEGLAB_plot.m) plots significant within-in and between-group differences by countouring significant pixels on the spectrogram. 



