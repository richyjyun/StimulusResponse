# StimulusResponse ([bioRxiv](https://www.biorxiv.org/content/10.1101/2022.03.30.486457v2))

Analysis of electrical stimulus (intracortical microstimulation) responses of single units. Data was collected using the Neural Interface Processor from Ripple Neuro with Utah arrays in the motor cortex of non-human primates. 

Example short-term excitatory response followed by long-term inhibitory response due to single-pulse intracortical microstimulation. 

<p align="center">
  <img width="700" height="500" src="https://github.com/richyyun/StimulusResponse/blob/main/InhibitionExample-01.png">
</p>

We tracked the changes in responses over time due to repetitive stimulation, their dependence on spontaneous firing rates and firing patterns, and the relationship between the excitatory and inhibitory responses.

## Analyses Performed
- Create a session list (SL) of all experimental sessions.
- Calculate both the excitatory and inhibitory responses for each detected spike. Excitatory response was calculated by thresholding the PSTH to find a peak. The inhibitory response was calculated by the time until the next spontaneous spike for each stimulus onset. 
- Compare whether pairs of spikes are more likely to be evoked by the same stimulus or not (cs_AnalyzePopulation). Also determine if each individual stimulus is more likely to evoke a population of spikes compared to random shuffling of spikes. 
- Determine changes in responses due to stimuls amplitude, distance from the stimulus site, and stimulus frequency.
- Determine whether the spontaneous firing rate is correlated with the excitatory or inhibitory responses over time
- Determine if the timing of each stimulus relative to the most recent spontanous spike affects the excitory or inhibitory response. 
