# VIP_reversal

Codes for analysis and making figures for the manuscript "Selective engagement of prefrontal VIP neurons in reversal learning".

## Authors
Jee Hyun Yi1,4, Young Ju Yoon1,4, Huijeong Jeong3, Seo Yeon Choe1, Min Whan Jung1,2,5*

1 Center for Synaptic Brain Dysfunctions, Institute for Basic Science, Daejeon 34141, Korea

2 Department of Biological Sciences, Korea Advanced Institute of Science and Technology, Daejeon 34141, Korea

3 Department of Neurology, University of California, San Francisco, CA 94158, USA.

4 These authors contributed equally to this work.

5 Lead contact *Correspondence: mwjung@kaist.ac.kr

## Code description
MATLAB version: R2017a

Each analysis takes a time scale of minutes.
*****
### MATLAB path setting
Make sure that <code>/functions</code> folder and data in <code>/data_VIP_CC</code> folder are on the MATLAB path.

No additional installations are required to run the code on the data.

*****
### Required MATLAB functions
Functions included in the <code>/functions</code> folder:
1. <code>othercolor</code> (Joshua Atkins (2011))
2. <code>cbrewer</code> (Charles Robert (2011))
3. <code>stateTime_zerofil</code>
4. <code>recorded_trial_types</code>
5. <code>behmat2eventlog</code>
6. <code>fun_ANCCR_3var_Ras</code>

ANCCR (https://github.com/namboodirilab/ANCCR) and CNMF-E (https://github.com/zhoupc/CNMF_E) packages are also required. (not included in this repository) (include them on your MATLAB path)

*****
### Behavioral and neural data
The data that support the findings of this study are available from the corresponding author upon reasonable request.


https://github.com/youngju0565/VIP_reversal
