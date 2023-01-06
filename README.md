# label_finder
in vivo isotope labling finder by comparing 12C vs 13C mass spectra

This program is designed to read two ascii files exported from QualBrowser, 12C and 13C data, and find 13C peaks in 13C data by matching to potentialmonoisotope peaks in 12C data. It is developed and tested for 13C in vivo isotope labeling of maize root tips. The citation will be made once the manuscript is accepted/published.

File list and explanation.
1. label_finder v7.py
: the main program. 

2. c12_1ave.txt
: averaged mass spectrum of longitudianl section of maize root tips grown in 12C-glucose

3. c13_1ave.txt
: averaged mass spectrum of longitudianl section of maize root tips grown in 13C-glucose

4. c12_1.txt
: averaged mass spectrum of cross-sectional section of maize root tips grown in 12C-glucose

5. c13_1.txt
: averaged mass spectrum of cross-sectional section of maize root tips grown in 13C-glucose
