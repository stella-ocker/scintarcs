# scintarcs

This repository contains example scripts used to analyze FAST filterbank data (pulsar tracking mode). These scripts were adapted for a specific data set and are not meant to be generally applicable (modifications will be required!). Contact me if you'd like to adapt these programs for your own usage cases. 

File descriptions:

foldpsrdat_3sec.py -- Program to fold filterbank data at given pulse period. Loads data, de-disperses to input DM, folds in 3 second subintegrations, downsamples in time to optimize storage. Saves 2D dynamic spectrum for bandpass calibration.

psrcal_tools.py -- Various utility functions for processing pulsar dynamic spectra.
