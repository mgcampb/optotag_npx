% Shuffle test to identify laser-responsive units
% MGC 11/10/2021

opt = struct;
opt.laser_offset = 2; % seconds relative to start of trial that laser pulses start
opt.ms_before = 10;
opt.ms_after = 40;
opt.bin_size_ms = 1;
opt.num_shuf = 1000;
opt.sigma_thresh = 5; % number of standard deviations above shuffle to be considered a significant bin

paths = struct;
paths.data = 'I:\My Drive\UchidaLab\DA_independence\neuropix_processed\ks3_thresh99';
paths.hgrk_analysis_tools = 'I:\My Drive\UchidaLab\code\HyungGoo';
addpath(genpath(paths.hgrk_analysis_tools));
paths.malcolm_functions = 'I:\My Drive\UchidaLab\code\malcolm_functions';
addpath(genpath(paths.malcolm_functions));
paths.save = 'I:\My Drive\UchidaLab\DA_independence\laser_response_shuffle_test';
if ~isfolder(paths.save)
    mkdir(paths.save);
end


%% get sessions
session = dir(fullfile(paths.data,'*.mat'));
session = {session.name}';
for i = 1:numel(session)
    session{i} = session{i}(1:end-4);
end

% subselect sessions to process
session = session(contains(session,'MC31'));


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

    laser_pulse_freq = dat.SessionData.TrialSettings(1).LaserPulseFrequency;
    laser_pulse_dur = dat.SessionData.TrialSettings(1).LaserPulseDuration;
    laser_pulse_num = dat.SessionData.TrialSettings(1).NumLaserPulse;

    num_pulse_total = nTrials * laser_pulse_num;

    laser_ts_all = sort(repmat(laser_ts,1,laser_pulse_num) + ...
        sort(repmat(1/laser_pulse_freq * (0:laser_pulse_num-1),1,numel(laser_ts))));

    %% Laser response versus shuffle (shifted spikes)

    fprintf('\tComputing PSTHs...\n');

    good_cells = dat.sp.cids(dat.sp.cgs==2);

    t = -opt.ms_before:opt.bin_size_ms:opt.ms_after;
    psth_laser = nan(numel(good_cells),nCond*laser_pulse_num,numel(t));
    psth_laser_shuf = nan(numel(good_cells),nCond*laser_pulse_num,numel(t),opt.num_shuf);

    tic
    for cIdx = 1:numel(good_cells)

        fprintf('\tCell %d/%d: cellID=%d\n',cIdx,numel(good_cells),good_cells(cIdx));
        spiket = 1000/opt.bin_size_ms * dat.sp.st(dat.sp.clu==good_cells(cIdx))';
        trigger = 1000/opt.bin_size_ms * laser_ts_all';
        grp = reshape(...
            repmat((TrialTypes-1)*laser_pulse_num,laser_pulse_num,1) + ...
            repmat((1:laser_pulse_num)',1,numel(TrialTypes)), ...
            numel(trigger),1);

        mint = 1000/opt.bin_size_ms * 20;
        maxt = 1000/opt.bin_size_ms * max(dat.sp.st);

        [~,psth] = plot_timecourse('timestamp',spiket,trigger,...
            trigger-opt.ms_before/opt.bin_size_ms,trigger+opt.ms_after/opt.bin_size_ms,grp,...
            'win_len',1,'resample_bin',1,'plot_type','none','grp_lim',nCond*laser_pulse_num);

        psth_laser(cIdx,:,:) = psth.mean;

        parfor shufIdx = 1:opt.num_shuf
            spiket_shuf = mod(spiket + unifrnd(mint,maxt),maxt);
            [~,psth] = plot_timecourse('timestamp',spiket_shuf,trigger,...
                trigger-opt.ms_before/opt.bin_size_ms,trigger+opt.ms_after/opt.bin_size_ms,grp,...
                'win_len',1,'resample_bin',1,'plot_type','none','grp_lim',nCond*laser_pulse_num);
            psth_laser_shuf(cIdx,:,:,shufIdx) = psth.mean;
        end

    end
    toc

    %% Identify significant cells

    mean_psth_laser_shuf = mean(psth_laser_shuf,4);
    std_psth_laser_shuf = std(psth_laser_shuf,[],4);
    thresh = mean_psth_laser_shuf + std_psth_laser_shuf * opt.sigma_thresh;
    psth_laser_thresholded = psth_laser>thresh;

    sig_laser_response = sum(psth_laser_thresholded(:,:,t>0 & t<=laser_pulse_dur*1000/opt.bin_size_ms),3)>0;

    %% Identify laser responsive cells

    % first pulse
    resp_first_pulse = nan(numel(good_cells),nCond);
    for cellIdx = 1:numel(good_cells)
        for condIdx = 1:nCond
            keep = 1 + (condIdx-1) * laser_pulse_num;
            resp_first_pulse(cellIdx,condIdx) = sig_laser_response(cellIdx,keep);
        end
    end

    % all pulses
    resp_all_pulse = nan(numel(good_cells),nCond);
    for cellIdx = 1:numel(good_cells)
        for condIdx = 1:nCond
            keep = (1:laser_pulse_num) + (condIdx-1) * laser_pulse_num;
            resp_all_pulse(cellIdx,condIdx) = all(sig_laser_response(cellIdx,keep));
        end
    end

    % at least one pulse
    resp_any_pulse = nan(numel(good_cells),nCond);
    for cellIdx = 1:numel(good_cells)
        for condIdx = 1:nCond
            keep = (1:laser_pulse_num) + (condIdx-1) * laser_pulse_num;
            resp_any_pulse(cellIdx,condIdx) = any(sig_laser_response(cellIdx,keep));
        end
    end

    %% save results
    save(fullfile(paths.save,opt.session),'opt','good_cells','psth_laser',...
        'mean_psth_laser_shuf','std_psth_laser_shuf','thresh',...
        'sig_laser_response','resp_first_pulse','resp_all_pulse','resp_any_pulse');
end
