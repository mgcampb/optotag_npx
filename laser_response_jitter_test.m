% Jitter test to identify laser-responsive units
% MGC 11/10/2021

% modified to jitter test MGC 11/10/2021

opt = struct;
opt.laser_offset = 2; % seconds relative to start of trial that laser pulses start
opt.ms_before = 10;
opt.ms_after = 40;
opt.num_jit = 500;
opt.jit_sd = 10; % in ms; jitter to add to each spike
opt.sigma_thresh = 5; % number of standard deviations above jitter to be considered a significant bin

paths = struct;
paths.data = 'I:\My Drive\UchidaLab\DA_independence\neuropix_processed\ks3_thresh99';
paths.libkm = 'C:\code\libkm';
addpath(genpath(paths.libkm));
paths.libmc = 'C:\code\libmc';
addpath(genpath(paths.libmc));
paths.save = 'I:\My Drive\UchidaLab\DA_independence\laser_response_jitter_test';
if ~isfolder(paths.save)
    mkdir(paths.save);
end


%% get sessions
session = get_mat_files(paths.data);

% subselect sessions to process
session = session(contains(session,'MC31'));

session = {'MC25_20211021'};


%% iterate over sessions
for sesh_idx = 1:numel(session)
    opt.session = session{sesh_idx};
    fprintf('Processing session %d/%d: %s\n',sesh_idx,numel(session),opt.session);

    dat = load(fullfile(paths.data,opt.session));

    %% get session data

    nTrials = dat.SessionData.nTrials;   
    TrialTypes = dat.SessionData.TrialTypes;
    laser_ts = dat.SessionData.TrialStartTimestamp + opt.laser_offset;
    if strcmp(dat.exp_params.bpod_protocol,'OdorLaserWater') % for this protocol, only keep rewarded trials (this includes rewarded non-laser trials as a control)
        keep = dat.SessionData.RewardedTrials == 1;
        TrialTypes = TrialTypes(keep);
        laser_ts = laser_ts(keep);
    end
    nCond = numel(unique(TrialTypes));

    PulseFreq = dat.SessionData.TrialSettings(1).LaserPulseFrequency;
    PulseDur = dat.SessionData.TrialSettings(1).LaserPulseDuration;
    NumPulse = dat.SessionData.TrialSettings(1).NumLaserPulse;

    num_pulse_total = nTrials * NumPulse;

    laser_ts_all = sort(repmat(laser_ts,1,NumPulse) + ...
        sort(repmat(1/PulseFreq * (0:NumPulse-1),1,numel(laser_ts))));

    %% Laser response versus jitter (shifted spikes)

    fprintf('\tComputing PSTHs...\n');

    good_cells = dat.sp.cids(dat.sp.cgs==2);

    t = -opt.ms_before:opt.ms_after;
    psth_laser = nan(numel(good_cells),nCond*NumPulse,numel(t));
    psth_laser_jit = nan(numel(good_cells),nCond*NumPulse,numel(t),opt.num_jit);

    tic
    for cIdx = 1:numel(good_cells)

        fprintf('\tCell %d/%d: cellID=%d\n',cIdx,numel(good_cells),good_cells(cIdx));
        spiket = 1000 * dat.sp.st(dat.sp.clu==good_cells(cIdx))';
        trigger = 1000 * laser_ts_all';
        grp = reshape(...
            repmat((TrialTypes-1)*NumPulse,NumPulse,1) + ...
            repmat((1:NumPulse)',1,numel(TrialTypes)), ...
            numel(trigger),1);

        [~,psth] = plot_timecourse('timestamp',spiket,trigger,...
            trigger-opt.ms_before,trigger+opt.ms_after,grp,...
            'win_len',1,'resample_bin',1,'plot_type','none','grp_lim',nCond*NumPulse);

        psth_laser(cIdx,:,:) = psth.mean;

        parfor jitIdx = 1:opt.num_jit
            spiket_jit = spiket + opt.jit_sd * randn(size(spiket));
            [~,psth] = plot_timecourse('timestamp',spiket_jit,trigger,...
                trigger-opt.ms_before,trigger+opt.ms_after,grp,...
                'win_len',1,'resample_bin',1,'plot_type','none','grp_lim',nCond*NumPulse);
            psth_laser_jit(cIdx,:,:,jitIdx) = psth.mean;
        end

    end
    toc

    %% Identify significant cells

    mean_psth_laser_jit = mean(psth_laser_jit,4);
    std_psth_laser_jit = std(psth_laser_jit,[],4);
    thresh = mean_psth_laser_jit + std_psth_laser_jit * opt.sigma_thresh;
    psth_laser_thresholded = psth_laser>thresh;

    sig_laser_response = sum(psth_laser_thresholded(:,:,t>0 & t<=PulseDur*1000),3)>0;

    %% Identify laser responsive cells

    % first pulse
    resp_first_pulse = nan(numel(good_cells),nCond);
    for cellIdx = 1:numel(good_cells)
        for condIdx = 1:nCond
            keep = 1 + (condIdx-1) * NumPulse;
            resp_first_pulse(cellIdx,condIdx) = sig_laser_response(cellIdx,keep);
        end
    end

    % all pulses
    resp_all_pulse = nan(numel(good_cells),nCond);
    for cellIdx = 1:numel(good_cells)
        for condIdx = 1:nCond
            keep = (1:NumPulse) + (condIdx-1) * NumPulse;
            resp_all_pulse(cellIdx,condIdx) = all(sig_laser_response(cellIdx,keep));
        end
    end

    % at least one pulse
    resp_any_pulse = nan(numel(good_cells),nCond);
    for cellIdx = 1:numel(good_cells)
        for condIdx = 1:nCond
            keep = (1:NumPulse) + (condIdx-1) * NumPulse;
            resp_any_pulse(cellIdx,condIdx) = any(sig_laser_response(cellIdx,keep));
        end
    end

    %% save results
    save(fullfile(paths.save,opt.session),'opt','good_cells','psth_laser',...
        'mean_psth_laser_jit','std_psth_laser_jit','thresh',...
        'sig_laser_response','resp_first_pulse','resp_all_pulse','resp_any_pulse');
end
