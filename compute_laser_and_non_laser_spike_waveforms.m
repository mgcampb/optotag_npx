% Computes average waveforms for laser-evoked and non-laser-evoked spikes
% using the SpikeGLX tool C_Waves
%
% MGC 10/12/2021

opt = struct;
opt.laser_offset = 2; % seconds relative to start of trial that laser pulses start
opt.ks_ver = 'ks3';

% cwaves options:
opt.samples_per_spike = 82;
opt.pre_samples = 30;
opt.num_spikes = 1000;
opt.snr_radius = 8;

paths = struct;
% data paths:
paths.data = 'F:\neuropix_processed\';
paths.raw_data_root = 'F:\ecephys\catGT\';
paths.dest_root = 'D:\waveforms\';
% code paths:
paths.hgrk_analysis_tools = 'C:\code\HyungGoo';
addpath(genpath(paths.hgrk_analysis_tools));
paths.malcolm_functions = 'C:\code\malcolm_functions';
addpath(genpath(paths.malcolm_functions));
paths.npy_matlab = 'C:\code\npy-matlab';
addpath(genpath(paths.npy_matlab));
paths.cwaves = 'C:\code\C_Waves-win\';

%% Get sessions
opt.session = dir(fullfile(paths.data,'*.mat'));
opt.session = {opt.session.name}';
for i = 1:numel(opt.session)
    opt.session{i} = opt.session{i}(1:end-4);
end

% subselect sessions to process
opt.session = opt.session(contains(opt.session,'MC30') | ...
    contains(opt.session,'MC31') | ...
    contains(opt.session,'MC33') | ...
    contains(opt.session,'MC34'));

%% Iterate over sessions
tic
for sesh_idx = 1:numel(opt.session)
    session_this = opt.session{sesh_idx};
    fprintf('\nAnalyzing session %d/%d: %s\n',sesh_idx,numel(opt.session),session_this);
    
    %% construct paths for this session
    paths.raw_data = [paths.raw_data_root 'catgt_' session_this '_g0\' session_this '_g0_imec0\'];
    paths.ks_output = [paths.raw_data 'imec0_' opt.ks_ver '\'];
    paths.spike_glx_bin = [paths.raw_data session_this '_g0_tcat.imec0.ap.bin'];
    paths.dest = [paths.dest_root session_this '\'];
    if ~isfolder(paths.dest)
        mkdir(paths.dest);
    end

    %% load data 
    fprintf('Loading data...\n');
    dat = load(fullfile(paths.data,session_this));
    good_cells = dat.sp.cids(dat.sp.cgs==2);
    cluster_info = tdfread(fullfile(paths.ks_output,'cluster_info.tsv'));

    %% get laser pulse intervals (all lasers)

    laser_pulse_duration = dat.SessionData.TrialSettings(1).LaserPulseDuration;
    laser_pulse_freq = dat.SessionData.TrialSettings(1).LaserPulseFrequency;
    num_laser_pulse = dat.SessionData.TrialSettings(1).NumLaserPulse;
    
    if strcmp(dat.exp_params.bpod_protocol,'OdorLaser')
        laser_trials = dat.SessionData.TrialTypes<=dat.SessionData.SettingsFile.NumLaser;
    elseif strcmp(dat.exp_params.bpod_protocol,'OdorLaserWater')
        laser_trials = dat.SessionData.TrialTypes<=dat.SessionData.SettingsFile.NumLaser & ...
            dat.SessionData.RewardedTrials==1;
    end
    laser_interval_base = dat.SessionData.TrialStartTimestamp(laser_trials) + opt.laser_offset;
    laser_interval_base = [laser_interval_base laser_interval_base+laser_pulse_duration];
    laser_bin_edges = [];
    for i = 1:num_laser_pulse
        laser_bin_edges = [laser_bin_edges laser_interval_base+(i-1)/laser_pulse_freq];
    end
    laser_bin_edges = sort(laser_bin_edges);

    %% find spikes within laser pulse (any laser) and spikes outside of laser pulse
    [~,~,laser_bin] = histcounts(dat.sp.st,laser_bin_edges); % make sure to use synchronized spikes (not uncorrected)
    laser_evoked = mod(laser_bin,2)==1;

    %% generate npy files - LASER_EVOKED

    fprintf('Generating npy files (laser-evoked spikes)...\n');

    paths.clus_table_npy = [paths.dest 'laser_evoked_clus_table.npy'];
    paths.clus_time_npy = [paths.dest 'laser_evoked_spike_times.npy'];
    paths.clus_lbl_npy = [paths.dest 'laser_evoked_spike_clusters.npy'];

    keep_spk = ismember(dat.sp.clu,good_cells) & laser_evoked;
    clus_time = uint64(dat.sp.st_uncorrected(keep_spk)*dat.sp.sample_rate); % make sure to use NON-synchronized spikes (uncorrected)
    clus_lbl = uint32(dat.sp.clu(keep_spk));

    % make clus_lbl correspond to rows of clus_table (ZERO BASED)
    for i = 1:numel(good_cells)
        clus_lbl(clus_lbl==good_cells(i)) = i-1;
    end

    clus_table = nan(numel(good_cells),2);
    for i = 1:numel(good_cells)
        clus_table(i,1) = sum(dat.sp.clu==good_cells(i) & laser_evoked);
        clus_table(i,2) = cluster_info.ch(cluster_info.id==good_cells(i));
    end
    clus_table = uint32(clus_table);

    writeNPY(clus_table,paths.clus_table_npy);
    writeNPY(clus_time,paths.clus_time_npy);
    writeNPY(clus_lbl,paths.clus_lbl_npy);

    %% generate and run cwaves command - LASER_EVOKED

    fprintf('Running C_Waves (laser-evoked spikes)...\n');

    command = sprintf(['%srunit.bat -spikeglx_bin=%s -clus_table_npy=%s -clus_time_npy=%s '...
        '-clus_lbl_npy=%s -dest=%s -samples_per_spike=%d -pre_samples=%d -num_spikes=%d -snr_radius=%d -prefix=laser_evoked'],...
        paths.cwaves, paths.spike_glx_bin, paths.clus_table_npy, paths.clus_time_npy, paths.clus_lbl_npy, paths.dest,...
        opt.samples_per_spike, opt.pre_samples, opt.num_spikes, opt.snr_radius);

    system(command);

    %% generate npy files - NON_LASER_EVOKED

    fprintf('Generating npy files (non-laser-evoked spikes)...\n');

    paths.clus_table_npy = [paths.dest 'not_laser_evoked_clus_table.npy'];
    paths.clus_time_npy = [paths.dest 'not_laser_evoked_spike_times.npy'];
    paths.clus_lbl_npy = [paths.dest 'not_laser_evoked_spike_clusters.npy'];

    keep_spk = ismember(dat.sp.clu,good_cells) & ~laser_evoked;
    clus_time = uint64(dat.sp.st_uncorrected(keep_spk)*dat.sp.sample_rate); % make sure to use NON-synchronized spikes (uncorrected)
    clus_lbl = uint32(dat.sp.clu(keep_spk));

    % make clus_lbl correspond to rows of clus_table (ZERO BASED)
    for i = 1:numel(good_cells)
        clus_lbl(clus_lbl==good_cells(i)) = i-1;
    end

    clus_table = nan(numel(good_cells),2);
    for i = 1:numel(good_cells)
        clus_table(i,1) = sum(dat.sp.clu==good_cells(i) & ~laser_evoked);
        clus_table(i,2) = cluster_info.ch(cluster_info.id==good_cells(i));
    end
    clus_table = uint32(clus_table);

    writeNPY(clus_table,paths.clus_table_npy);
    writeNPY(clus_time,paths.clus_time_npy);
    writeNPY(clus_lbl,paths.clus_lbl_npy);

    %% generate and run cwaves command - NON_LASER_EVOKED

    fprintf('Running C_Waves (non-laser-evoked spikes)...\n');

    command = sprintf(['%srunit.bat -spikeglx_bin=%s -clus_table_npy=%s -clus_time_npy=%s '...
        '-clus_lbl_npy=%s -dest=%s -samples_per_spike=%d -pre_samples=%d -num_spikes=%d -snr_radius=%d -prefix=not_laser_evoked'],...
        paths.cwaves, paths.spike_glx_bin, paths.clus_table_npy, paths.clus_time_npy, paths.clus_lbl_npy, paths.dest,...
        opt.samples_per_spike, opt.pre_samples, opt.num_spikes, opt.snr_radius);

    system(command);
    toc
end