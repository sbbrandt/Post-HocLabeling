function [epo_sources, Ax] = extract_sources(epo_narrow,type, options)
switch type
    case 'hd'
        [epo_sources, Ax] = extract_sourceHD(epo_narrow, options);
    case 'ica'
        [epo_sources, Ax] = extract_sourceICA(epo_narrow, options);
end
end

function [epo_ica_mara, Ax] = extract_sourceICA(epo_narrow, options)

cnt_narrow = epo_narrow;
cnt_narrow.x = custom_epo2cnt(epo_narrow.x,'Ne',size(epo_narrow.x,3));
[cnt_ica, W_ica, A_ica] = proc_fastICA(cnt_narrow, 'fasticaParams',...
    {'verbose', 'on','displayMode','off','numOfIC',options.N_compICA});

epo_ica = cnt_ica;
epo_ica.x = custom_epo2cnt(cnt_ica.x,'Ne',size(epo_narrow.x,3));

[~, info] = proc_MARA(cnt_ica,epo_ica.clab,A_ica);
goodcomp = find(info.out_p < 1e-8);
epo_ica_mara = proc_selectChannels(epo_ica, goodcomp);

Ax = struct;
Ax.Ax_all = A_ica(:,goodcomp);
Ax.W_ica = W_ica(:,goodcomp);
Ax.clab = epo_narrow.clab;
end

function [epo_sources, Ax] = extract_sourceHD(epo_narrow, options)
[Nc, Nt, Ne] = size(epo_narrow.x);

x_concat = custom_epo2cnt(epo_narrow.x);

[epo_narrow, Ax, inv_operator] = ...
            load_sr_model(epo_narrow,  options);

x_concat = custom_epo2cnt(epo_narrow.x);

s_allNoEnv = single(x_concat)*single(inv_operator)';
s_allNoEnv = custom_epo2cnt(s_allNoEnv,'Ne',Ne);
epo_sources = epo_narrow;
epo_sources.x = s_allNoEnv;
epo_sources.clab = strread(num2str(1:size(s_allNoEnv,2)),'%s');

end

function [bbci_eeg, Ax, inv_operator, cortex, fm] = load_sr_model(bbci_eeg,  options)
load(options.file_headModel);

% check channels contained in head model and EEG
bbci_eeg = proc_selectChannels(bbci_eeg,Ax.clab);
idx_matches = find(ismember(Ax.clab,bbci_eeg.clab));
Ax.Ax_all = Ax.Ax_all(idx_matches,:);
Ax.clab = Ax.clab(idx_matches);
Ax.mnt = mnt_setElectrodePositions(Ax.clab);

L = Ax.Ax_all;
[Nc,Ns] = size(L);

%optimize regularization parameter
x_cnt = custom_epo2cnt(bbci_eeg.x);
x_cnt = x_cnt(randsample(size(x_cnt,1),0.1*size(x_cnt,1)),:);
opt_fun = @(lambda) gcv(lambda,L,x_cnt',speye(Ns,Ns),speye(Nc,Nc));
reg_par = fminbnd(opt_fun,0,1);

% Generate inverse operator for source reconstruction
% reg_par = 1e-4;
Qe = reg_par*eye(Nc,Nc);
Q = speye(Ns,Ns);
inv_operator = Q*L'/(Qe + L*Q*L');


end

function val = gcv(reg_par,L,Y,Q,I_Nc)
Qe = reg_par*I_Nc;
inv_operator = Q*L'/(Qe + L*Q*L');
X = inv_operator*Y;
num = norm(L*X -Y)^2;
den = trace(I_Nc-L*inv_operator)^2;
val = num/den;
end
