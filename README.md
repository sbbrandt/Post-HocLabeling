# Post-hoc labeling of EEG for benchmarking of decoding methods
This is the source code for the post-hoc labelled dataset generation framework (preprint available at http://arxiv.org/abs/1711.08208)


## Getting started
1. Run the download_dependecies.m script to download and configure the bbci_toolbox.
2. Download EEG data available at https://zenodo.org/record/1065107#.WhaCX3XyvCI and place the .mat files under ./data
3. Run the script test_postHocLabelling.m to see usage example.

## Comments on the New York's head model -
The Original head model can be downloaded at https://www.parralab.org/nyhead/
The .mat available here was stripped to the strictly necessary data in order to run this toolbox.
