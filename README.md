# msa-online

Online *microstate sequence analysis (msa)* in 3 steps:

1. Data input
2. Parse data
3. Analyze

**Data input**  
- copy & paste microstate sequence as text
- load microstate sequence from file

**Parse**  
Microstate symbols will be found automatically, supported symbols are lower case letters (a-z), upper case letters (A-Z) and integers (0-9).  

**Analyze**  
* Distribution of symbols
* Transition probability matrix
* Zero-order Markov test
* First-order Markov test
* Shannon entropy of the sequence
* Finite entropy rate
* Active information storage
* Partial autoinformation function
* Autoinformation function

These analyses were used in the publications:

[1] F. von Wegner, S. Bauer, F. Rosenow, J. Triesch, H. Laufs, “EEG microstate periodicity explained by rotating phase patterns 
    of resting-state alpha oscillations.”, NeuroImage, 224:117372, 2021.  
[2] F. von Wegner, H. Laufs, “Information-theoretical analysis of EEG microstate sequences in Python.”, Front Neuroinform, 12:30, 2018.  
[3] F. von Wegner, P. Knaut, H. Laufs, “EEG microstate sequences from different clustering algorithms are information-theoretically 
    invariant.”, Front Comp Neurosci, 12:70, 2018.  
[4] F. von Wegner, E. Tagliazucchi, H. Laufs, “Information theoretical analysis of resting state EEG microstate sequences - 
    non-Markovianity, non-stationarity and periodicities.”, NeuroImage, 158:99-111, 2017.
