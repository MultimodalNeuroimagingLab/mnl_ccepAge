# mnl_ccepAge

This repository contains the scripts and functions necessary to reproduce the analyses and figures in this manuscript:

- The development of efficient communication in the human connectome. D. van Blooijs, J.F. van der Aar, G.J.M. Huiskamp, G. Castegnaro, M. Demuru, W.J.E.M. Zweiphenning, P. van Eijsden, K. J. Miller, F.S.S. Leijten, D. Hermes bioRxiv 2022.03.14.484297; doi: https://doi.org/10.1101/2022.03.14.484297


# Generating the figures
Scripts to process the data, detect N1 responses and add tract information:
- scripts/ccep01_averageCCEPs.m
- scripts/ccep02_loadN1.m
- scripts/ccep03_addtracts.m

Scripts to make the figure panels:
- makeFig1A_plotMNI.m
- makeFig1B_ExampleResponse.m
- makeFig2_plotResponsesAge.m
- makeFig3and4_plotN1_meanacrossage.m


To make all functions work, an m-file called personalDataPath.m should be stored in the root dir. This file should have the following content:
```
function localDataPath = personalDataPath()
% function that contains local data path, is ignored in .gitignore
localDataPath.input = '/my/path/to/load/data/';
localDataPath.output = '/my/path/to/save/data/';
addpath('/my/path/to/fieldtrip')
addpath('/my/path/to/leadDBS')
ft_defaults
```

# Dependencies
- Fieldtrip: http://www.fieldtriptoolbox.org/  
- Lead DBS: https://www.lead-dbs.org/


# Acknowledgements
The project is funded by the National Institute Of Mental Health of the National Institutes of Health under Award Number R01MH122258 to Dora Hermes (Mayo Clinic) and by the Dutch Epilepsy Foundation (NEF #17-07).


# Contact
----------------------------
Please contact Dora Hermes (hermes.dora@mayo.edu), Dorien van Blooijs or Max van den Boom for questions. 
