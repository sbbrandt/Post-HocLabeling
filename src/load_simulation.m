function [epo_sources,Ax] = load_simulation(eeg_data, varargin)
% load_simulation - Post-hoc labelled dataset generation
%
% Generates target variables of hidden oscillatory sources from arbitrary
% EEG recordings
%
%Synopsis:
% [epo_sources,Ax] = load_simulation(eeg_data, <OPT>)
%
%Arguments:
% eeg_data    - data structure containing:
%               cnt: Data structure with continuous EEG data (BBCI_toolbox
%               format)
%               vmrk: (Arbitrary) markers to segment cnt data (BBCI_toolbox
%               format)
%               iart: Indices of markers in vmrk which are to be rejected
%               due to artifacts
%
% OPT - struct or property/value list of optional properties:
%  .file_headModel - full file name of the head model file
%  .ford - order of the butterworth filter to perform frequency filtering
%  .cFreq - central frequency of the target oscillatory source
%  .wFreq - width of the band of the target oscillatory source
%  .winddowLength - length of the segments to extract from cnt
%  .type - Type of label generation (can be 'hd' for head model or 'ica')
%  .select_sources - number of source to return, 'all' for return all
%  sources computed

%Returns:
% epo_sources  - segmented version of the envelop of the sources, to be
% used as labels
% Ax           - spatial patterns corresponding to each of the sources in
% epo_sources
% 
% References:
%
% soon

global BTB

if ~isstruct(varargin)
    options = propertylist2struct(varargin{:});
end

options = set_defaults(options,...
    'file_headModel', fullfile('data/sa_nyhead_simplified.mat'), ... % location of the head 
    'ford', 5, ...
    'cFreq', [10] , ...
    'wFreq', [4], ...
    'windowLength', 1000, ...
    'N_compICA', 20,...
    'type', 'hd', ...
    'select_sources', 'all');

cnt = eeg_data.cnt;
iart = eeg_data.iart;
vmrk = eeg_data.vmrk;
clear eeg_data

%% Generate filters for target frequency bands and filter EEG 
fband = [options.cFreq-options.wFreq/2, ...
        options.cFreq+options.wFreq/2];
[zfilt,pfilt,kfilt] = butter(options.ford,...
                fband/(cnt.fs/2));

cnt_narrow = proc_filt(cnt,zfilt,pfilt,kfilt);

%% Segment continuous data
vmrk_org = vmrk;
vmrk = mrk_selectEvents(vmrk,'not',iart);
epo_narrow = proc_segmentation(cnt_narrow, vmrk, ....
            [0 options.windowLength]);

%% Compute sources
[epo_sources, Ax] = extract_sources(epo_narrow, options.type, options);
Ns = size(epo_sources.x, 2);

%% Select sources (or return all)
if ischar(options.select_sources)
    idx_tSource = 1:Ns;
elseif isnumeric(options.select_sources)
    idx_tSource = options.select_sources;    
    epo_sources = proc_selectEpochs(epo_sources, idx_tSource);
end

%% Obtain envelope
Ne = size(epo_sources.x,3);
for idx_s = 1:length(idx_tSource)
     cs_target = custom_epo2cnt(epo_sources.x(:,idx_s,:));
     cs_target = abs(hilbert(cs_target));
     epo_sources.x(:,idx_s,:) = custom_epo2cnt(cs_target,'Ne',Ne);
end     



end
