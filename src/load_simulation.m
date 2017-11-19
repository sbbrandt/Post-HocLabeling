function [epo_sources,Ax] = load_simulation(eeg_data, varargin)

global BTB

if ~isstruct(varargin)
    options = propertylist2struct(varargin{:});
end

options = set_defaults(options,...
    'file_headModel', fullfile('data/sa_nyhead_simplified.mat'), ...
    'ford', 5, ...
    'cFreq', [10] , ...
    'wFreq', [4], ...
    'windowLength', 1000, ...
    'eventSpacing', 1000,...
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
