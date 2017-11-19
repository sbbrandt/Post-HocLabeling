% Sample script for generation of labeled dataset from pre-recorded
% EEG signals.
%
% sebastian.castano@blbt.uni-freiburg.de
% 17th November 2017

clear; close all; clc
addpath('src')

%% load eeg data and generate set of possible target sources (post-hoc labeling)
subject_str = 'S7';
options = struct(...
    'file_headModel', fullfile('data/sa_nyhead_simplified.mat'), ...
    'ford', 5, ...
    'cFreq', [10] , ...
    'wFreq', [4], ...
    'windowLength', 1000, ...
    'N_compICA', 20,...
    'type', 'ica', ...
    'select_sources', 'all');

eeg_data = load(fullfile('data',subject_str));
[epo_sources, Ax] = load_simulation(eeg_data, options);

%% filter EEG data to target frequency bands
[zfilt,pfilt,kfilt] = butter(5,[8,12]/(eeg_data.cnt.fs/2));
cnt_filt = proc_filt(eeg_data.cnt,zfilt,pfilt,kfilt);
epo = proc_segmentation(cnt_filt, eeg_data.vmrk, [0,1000]);
epo = proc_selectEpochs(epo, 'not', eeg_data.iart);

%% Train and apply SPoC
close all
% select random source as target and extract labels
ix_targetIndex = randi(size(epo_sources.x,2));
epo_target = proc_selectChannels(epo_sources,ix_targetIndex);
z = squeeze(mean(epo_target.x,1));
epo.y = z';

% split train/test set
Ne = size(epo.x,3);
[ix_train,ix_val,~] = divideblock(Ne,0.7,0.3,0);
epo_tr = proc_selectEpochs(epo, ix_train, 'RemoveVoidClasses', 0);
epo_val = proc_selectEpochs(epo, ix_val,  'RemoveVoidClasses', 0);

% train SPoC
[~, W, A_est] = proc_spoc(epo_tr);
w = W(:,1);
a_zest = A_est(:,1);

% apply SPoC
epo_targetPred= proc_linearDerivation(epo_val, W(:,1), 'prependix','spoc');

% get performance measurements

% correlation rho
z_val = z(ix_val);
z_pred = squeeze(var(epo_targetPred.x,[],1));
rho = corrcoef(z_val,z_pred);
rho = rho(1,2);
fprintf('Performance: Correlation rho=%.4f\n ',rho)

% angle between spatial patterns (first order channels so that both spatial
% patterns are comparable)
[~,ix_order_channels] = ismember(epo.clab,Ax.clab);
val_ch = find(ix_order_channels ~= 0);
ix_order_channels = ix_order_channels(val_ch);

a_z = Ax.Ax_all(ix_order_channels,ix_targetIndex);
alpha = acos(dot(a_zest(val_ch),a_z)/(norm(a_z)*norm(a_zest)));
if (alpha > pi/2)
    alpha = pi- alpha;
end
fprintf('Performance: Angle alpha=%.4f\n ',alpha)

% plot original and estimated spatial patterns
mnt = mnt_setElectrodePositions(Ax.clab(ix_order_channels));
figure;
subplot(1,2,1); plot_scalp(mnt, a_z, 'ScalePos', 'none'); 
title('original sp')

subplot(1,2,2); plot_scalp(mnt, double(a_zest(val_ch)), 'ScalePos', 'none'); 
title('estimated sp')



