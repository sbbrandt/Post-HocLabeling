function auxiliary_generate4Python(cfg_file)
addpath('~/libraries')
set_localpaths();
cfg = ini2struct(cfg_file);
if strcmp(fieldnames(cfg),'config')
    cfg = cfg.config;
end
%make str entries num
fn = fieldnames(cfg);
for ii=1:numel(fn)
    if ~isempty(str2num(cfg.(fn{ii})))
        cfg.(fn{ii}) = str2num(cfg.(fn{ii}));
    end
end

assert(diff(cfg.fband) > 0)
fc = mean(cfg.fband);
fw = diff(cfg.fband);
winlength = cfg.winlength;
type = cfg.type;
power_quantile = cfg.power_quantile;
discrete_quantile = cfg.discrete_quantiles;
eeg_data_file=cfg.eeg_data_file;
file_name_save=cfg.file_name_save;
%Ne = cfg.ne;
%eta = cfg.eta;

%% get random from .../S1-S7, if path is .../random
%if ismember('random',split(eeg_data_file,'/'))
%    eeg_data_file = join([eeg_data_file(1:end-6),'S',string(randi(7))],'');
%end

%%
options = struct(...
    'file_headModel', fullfile('data/sa_nyhead_simplified.mat'), ...
    'ford', 5, ...
    'cFreq', fc , ...
    'wFreq', fw, ...
    'windowLength', winlength, ...
    'N_compICA', 20,...
    'type', type, ...
    'select_sources', 'all', ...
    'power_quartile',power_quantile,...
    'discrete_quantiles', discrete_quantile);
    %'select_sources', 'all');

window_offset = [0,0];

eeg_data = load(eeg_data_file);

% skip the first and last n events
nskip = 0;
ntotal = length(eeg_data.vmrk.y);
%eeg_data.vmrk = mrk_selectEvents(eeg_data.vmrk,'not',[1:nskip, ntotal-nskip:ntotal]);
%eeg_data.iart = eeg_data.iart-nskip;

[epo_sources, Ax] = load_simulation(eeg_data, options);

%%
vmrk = mrk_selectEvents(eeg_data.vmrk, 'not', eeg_data.iart);
[epo, complete]= proc_segmentation(eeg_data.cnt, vmrk, [window_offset(1),winlength + window_offset(2)]);

epo_sources = proc_selectEpochs(epo_sources,complete);

cnt=eeg_data.cnt;
vmrk = mrk_selectEvents(vmrk,complete);
% select data length
% if Ne~=0
%     epo_sources.x = epo_sources.x(:,:,1:Ne);
%     epo_sources.y = epo_sources.y(1:Ne);
%     epo_sources.event.desc(end+1-Ne:end);
%     epo_sources.mrk_info.event.desc(end+1-Ne:end);
%     vmrk.y = vmrk.y(1:Ne);
%     vmrk.time = vmrk.time(1:Ne);
%     vmrk.event.desc(end+1-Ne:end);
% end
% add noise
% if eta~=0
%     epo_sources.x = epo_sources.x + (1-(1-eta)^2)/(1-eta)^2 ...
%         *var(mean(epo_sources.x))*randn(1,1,size(epo_sources.x,3));
% end
%%
save(file_name_save,'cnt','vmrk','epo_sources','window_offset', 'winlength');
%ntrials = length(eeg_data.vmrk.time);
%eeg_data.iart = unique([1,2,n_trials-1,n_trials-2,eeg_data.iart])

%[epo, complete]= proc_segmentation(eeg_data.cnt, vmrk, [-500,winlength+500]);
return