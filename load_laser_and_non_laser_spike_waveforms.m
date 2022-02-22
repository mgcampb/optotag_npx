% loads the results of compute_laser_and_non_laser_spike_waveforms.m
% and saves them to the processed data files
% MGC 11/10/2021

opt = struct;

paths = struct;
paths.data = 'F:\neuropix_processed';
paths.cwaves_output = 'D:\waveforms\';
paths.npy_matlab = 'C:\code\npy-matlab';
addpath(genpath(paths.npy_matlab));
paths.save = 'D:\waveforms\waveform_correlation\';


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

    dat = load(fullfile(paths.data,opt.session));
    good_cells = dat.sp.cids(dat.sp.cgs==2);
    
    laser_evoked_waveforms_all_chan = readNPY(fullfile(paths.cwaves_output,opt.session,'laser_evoked_mean_waveforms.npy'));
    non_laser_evoked_waveforms_all_chan = readNPY(fullfile(paths.cwaves_output,opt.session,'not_laser_evoked_mean_waveforms.npy'));

    [~,maxidx] = max(range(non_laser_evoked_waveforms_all_chan,3),[],2);
    
    laser_evoked_waveforms = nan(numel(good_cells),size(laser_evoked_waveforms_all_chan,3));
    non_laser_evoked_waveforms = nan(numel(good_cells),size(laser_evoked_waveforms_all_chan,3));
    for i = 1:numel(good_cells)
        laser_evoked_waveforms(i,:) = laser_evoked_waveforms_all_chan(i,maxidx(i),:);
        non_laser_evoked_waveforms(i,:) = non_laser_evoked_waveforms_all_chan(i,maxidx(i),:);
    end
    
    waveform_corr = diag(corr(laser_evoked_waveforms',non_laser_evoked_waveforms'));

    % save data
    save(fullfile(paths.save,opt.session),'opt','good_cells',...
        'laser_evoked_waveforms','non_laser_evoked_waveforms','waveform_corr');
end