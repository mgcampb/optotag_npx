% load qc metrics and waveforms from file
% MGC 11/10/2021

opt = struct;

paths = struct;
paths.data = 'F:\neuropix_processed'; % location of processed npx data files
paths.catGT = 'F:\ecephys\catGT\';
paths.ks_ver = 'ks3';
paths.npy_matlab = 'C:\code\npy-matlab';
addpath(genpath(paths.npy_matlab));


%% get sessions
session = dir(fullfile(paths.data,'*.mat'));
session = {session.name}';
for i = 1:numel(session)
    session{i} = session{i}(1:end-4);
end

% subselect sessions to process
session = session(contains(session,'MC30') | ...
    contains(session,'MC31') | ...
    contains(session,'MC33') | ...
    contains(session,'MC34'));

%% iterate over sessions
for sesh_idx = 1:numel(session)
    opt.session = session{sesh_idx};
    fprintf('Processing session %d/%d: %s\n',sesh_idx,numel(session),opt.session);
    paths.ks_dir = [paths.catGT 'catgt_' opt.session '_g0\' opt.session '_g0_imec0\imec0_' paths.ks_ver];
    
    % load spike data 
    dat = load(fullfile(paths.data,opt.session));
    assert(issorted(dat.sp.cids));

    % load metrics (qc and waveform)
    metrics_file = dir(fullfile(paths.ks_dir,'metrics*.csv'));
    metrics_file = sort({metrics_file.name}');
    metrics_file = metrics_file{end};
    
    metrics = readtable(fullfile(paths.ks_dir,metrics_file));
    metrics = metrics(ismember(metrics.cluster_id,dat.sp.cids),:);
    assert(all(metrics.cluster_id==dat.sp.cids'));

    % load mean waveforms
    waveforms_file = dir(fullfile(paths.ks_dir,'mean_waveforms*.npy'));
    waveforms_file = sort({waveforms_file.name}');
    waveforms_file = waveforms_file{end};

    mean_waveforms = readNPY(fullfile(paths.ks_dir,waveforms_file));
    mean_waveforms = mean_waveforms(dat.sp.cids+1,:,:); % cluster IDs start at 0 (hence the + 1)

    % replace fields
    if isfield(dat.sp,'waveform_metrics')
        dat.sp = rmfield(dat.sp,'waveform_metrics');
    end
    if isfield(dat.sp,'qc_metrics')
        dat.sp = rmfield(dat.sp,'qc_metrics');
    end
    if isfield(dat.sp,'mean_waveforms_cluster_id')
        dat.sp = rmfield(dat.sp,'mean_waveforms_cluster_id');
    end

    dat.sp.mean_waveforms = mean_waveforms;
    dat.sp.metrics = metrics;

    % save updated data
    sp = dat.sp;
    save(fullfile(paths.data,opt.session),'sp','-append');
end