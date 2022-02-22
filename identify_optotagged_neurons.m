% identify optotagged neurons using:
%   1) QC metrics to filter out low-quality units
%   2) shuffle test to identify laser-responsive units
%   3) waveform correlation (laser- vs non-laser-evoked spikes)
%   4) collision test to confirm antidromic stimulation
% MGC 11/10/2021

% generates "cell_table" containing all this info and more
% MGC 11/11/2021

opt = struct;

opt.make_figs = false;

opt.isi_viol_thresh = 0.5;
opt.amplitude_cutoff_thresh = 0.1;
opt.presence_ratio_thresh = 0.9;
opt.waveform_corr_thresh = 0.9;
opt.min_num_resp_pulse = 5;

opt.laser_offset = 2; % seconds relative to start of trial that laser pulses start
opt.ms_before = 100; % for laser psth
opt.ms_after = 600;
opt.bin_size_ms = 1;

opt.peak_fr_cutoff = 15; % max fr in first pulse (when deemed laser-responsive) must be bigger than this
opt.peak_fr_zscore_cutoff = 5; % max fr in first pulse (when deemed laser-responsive) must exceed baseline (10 ms before laser pulse) by this many z-scores

% collision test params
opt.collision_test_window = [-7 -2]; % in ms, window before first significant bin in which to look for trials with spontaneous spikes

paths = struct;
paths.data = 'F:\neuropix_processed\';
paths.wv_corr = 'D:\waveforms\waveform_correlation\';
paths.shuf_test = 'D:\laser_response_shuffle_test\';
paths.hgrk_analysis_tools = 'C:\code\HGRK_analysis_tools\';
addpath(genpath(paths.hgrk_analysis_tools));
% paths.figs = [paths.gdrive 'figs\identify_optotagged_neurons\optotagged_only'];
paths.figs = 'C:\figs\da_independence\identify_optotagged_neurons\optotagged_only\';
if opt.make_figs && ~isfolder(paths.figs)
    mkdir(paths.figs);
end


%% get sessions
session = dir(fullfile(paths.data,'*.mat'));
session = {session.name}';
for i = 1:numel(session)
    session{i} = session{i}(1:end-4);
end

session = session(contains(session,'MC28'));

%% gather cells from all sessions

max_size = 10000;
cell_table = table();
cell_table.Mouse = cell(max_size,1);
cell_table.Session = cell(max_size,1);
cell_table.CellID = nan(max_size,1);
cell_table.UniqueID = cell(max_size,1);
cell_table.mean_fr = nan(max_size,1);
cell_table.isi_viol = nan(max_size,1);
cell_table.presence_ratio = nan(max_size,1);
cell_table.amplitude_cutoff = nan(max_size,1);
cell_table.waveform = nan(max_size,82);
cell_table.peak_channel = nan(max_size,1);
cell_table.waveform_amplitude = nan(max_size,1);
cell_table.waveform_corr = nan(max_size,1);
cell_table.ShufTest = nan(max_size,40);
cell_table.RespAllPulse = nan(max_size,4);
cell_table.RespAnyPulse = nan(max_size,4);
cell_table.RespFirstPulse = nan(max_size,4);
cell_table.LaserResponse = nan(max_size,4);
cell_table.PassFRCutoff = nan(max_size,4);
cell_table.FirstSigBin = nan(max_size,40);
cell_table.PeakResponseBin = nan(max_size,40);
cell_table.CollisionTest = nan(max_size,1); 
cell_table.Optotagged = nan(max_size,1);

ctr = 1;
tic;
for sesh_idx = 1:numel(session)
    opt.session = session{sesh_idx};
    fprintf('Processing session %d/%d: %s\n',sesh_idx,numel(session),opt.session);

    % load data
    dat = load(fullfile(paths.data,opt.session));
    good_cells = dat.sp.cids(dat.sp.cgs==2);
    nCells = numel(good_cells);
    TrialTypes = dat.SessionData.TrialTypes;
    nCond = numel(unique(TrialTypes));
    nPulse = dat.SessionData.TrialSettings(1).NumLaserPulse;
    pulseDurMs = dat.SessionData.TrialSettings(1).LaserPulseDuration * 1000;
    laser_freq = dat.SessionData.TrialSettings(1).LaserPulseFrequency;
    laser_ts = dat.SessionData.TrialStartTimestamp + opt.laser_offset;    
    laser_ts_all = sort(repmat(laser_ts,1,nPulse) + ...
        sort(repmat(1/laser_freq * (0:nPulse-1),1,numel(laser_ts))));
    
    if strcmp(dat.exp_params.bpod_protocol,'OdorLaser')
        laser_target = [dat.exp_params.laser_target {'none'}];
    elseif strcmp(dat.exp_params.bpod_protocol,'OdorLaserWater')
        laser_target = [dat.exp_params.laser_target {'none'} {'none'}];
    end
    
    % reorder laser response data to be {VS, DMS, DLS, none}
    laser_reorder = (1:nCond)';
    laser_reorder(1) = find(strcmp(laser_target,'VS'));
    laser_reorder(2) = find(strcmp(laser_target,'DMS'));
    laser_reorder(3) = find(strcmp(laser_target,'DLS'));
    
    % reorder individual laser pulses
    pulse_reorder = reshape((nPulse*repmat(laser_reorder-1,1,nPulse)+repmat(1:nPulse,nCond,1))',nPulse*nCond,1);
    
    % reorder TrialTypes
    TrialTypes2 = TrialTypes;
    for i = 1:nCond
        TrialTypes2(TrialTypes==laser_reorder(i)) = i;
    end
    TrialTypes = TrialTypes2;
    
    % enter cell IDs into table
    idx = ctr:ctr+nCells-1;
    ctr = ctr+nCells;
    mouse_this = opt.session(1:4);
    cell_table.Mouse(idx) = {mouse_this};
    cell_table.Session(idx) = {opt.session};
    cell_table.CellID(idx) = good_cells;
    for i = 1:nCells
        cell_table.UniqueID{idx(i)} = sprintf('%s_c%d',opt.session,good_cells(i));
    end
    
    % cell metrics
    keep_cell = dat.sp.cgs==2;
    metrics = dat.sp.metrics(keep_cell,:);
    cell_table.mean_fr(idx) = metrics.firing_rate;
    cell_table.isi_viol(idx) = metrics.isi_viol;
    cell_table.presence_ratio(idx) = metrics.presence_ratio;
    cell_table.amplitude_cutoff(idx) = metrics.amplitude_cutoff;
    cell_table.peak_channel(idx) = metrics.peak_channel;
    cell_table.waveform_amplitude(idx) = metrics.amplitude;
    
    % avg waveform
    mean_waveforms = dat.sp.mean_waveforms(keep_cell,:,:);
    [~,maxidx] = max(range(mean_waveforms,3),[],2);
    wv = nan(nCells,82);
    for i = 1:nCells
        wv(i,:) = mean_waveforms(i,maxidx(i),:);
    end
    cell_table.waveform(idx,:) = wv;
    
    % load waveform correlation (laser vs non-laser)
    wv_corr = load(fullfile(paths.wv_corr,opt.session));
    cell_table.waveform_corr(idx) = wv_corr.waveform_corr;
    
    % load shuffle test results
    shuf_test = load(fullfile(paths.shuf_test,opt.session));
    cell_table.ShufTest(idx,:) = shuf_test.sig_laser_response(:,pulse_reorder);
    cell_table.RespAllPulse(idx,:) = shuf_test.resp_all_pulse(:,laser_reorder);
    cell_table.RespAnyPulse(idx,:) = shuf_test.resp_any_pulse(:,laser_reorder);
    cell_table.RespFirstPulse(idx,:) = shuf_test.resp_first_pulse(:,laser_reorder);

    
    % laser response (shuffle test)
    % laser_responsive = any(shuf_test.resp_all_pulse,2);
    num_laser_resp_pulse = nan(numel(good_cells),nCond);
    for i = 1:nCond
        keep_idx = (1:nPulse) + (i-1) * nPulse;
        num_laser_resp_pulse(:,i) = sum(shuf_test.sig_laser_response(:,keep_idx),2);
    end
    laser_response = shuf_test.resp_first_pulse & ...
        num_laser_resp_pulse>=opt.min_num_resp_pulse;
    laser_response = laser_response(:,laser_reorder);
    cell_table.LaserResponse(idx,:) = laser_response;
    
    % enforce some additional FR cutoffs (there were lots of low-FR cells
    % that passed the shuffle test):
    t_psth = -shuf_test.opt.ms_before:shuf_test.opt.bin_size_ms:shuf_test.opt.ms_after;
    baseline_fr = mean(shuf_test.psth_laser(:,1:nPulse:end,t_psth<0),3);
    baseline_std = std(shuf_test.psth_laser(:,1:nPulse:end,t_psth<0),[],3);
    max_fr = max(shuf_test.psth_laser(:,1:nPulse:end,t_psth>0 & t_psth<=pulseDurMs),[],3);
    pass_fr_cutoff = max_fr > opt.peak_fr_cutoff & ...
        max_fr > baseline_fr + baseline_std * opt.peak_fr_zscore_cutoff;
    pass_fr_cutoff = pass_fr_cutoff(:,laser_reorder);
    cell_table.PassFRCutoff(idx,:) = pass_fr_cutoff;

    
    % collision test
    collision_test = nan(nCells,1);
    optotagged = nan(nCells,1);
    for cIdx = 1:nCells
        fprintf('\tCell %d/%d: cellID=%d\n',cIdx,nCells,good_cells(cIdx));
        
        spiket = 1000/opt.bin_size_ms * dat.sp.st(dat.sp.clu==good_cells(cIdx))';
        trigger = 1000/opt.bin_size_ms * laser_ts';
        
        sig_resp_laser_this = squeeze(shuf_test.psth_laser(cIdx,:,t_psth>0 & t_psth<=pulseDurMs) > ...
            shuf_test.thresh(cIdx,:,t_psth>0 & t_psth<=pulseDurMs));
        sig_resp_laser_this = sig_resp_laser_this(pulse_reorder,:);

        first_sig_bin = nan(nCond * nPulse,1);

        for condIdx = 1:numel(first_sig_bin)
            first_bin_this = find(sig_resp_laser_this(condIdx,:)>0,1,'first');
            if ~isempty(first_bin_this)
                first_sig_bin(condIdx) = first_bin_this;
            end
        end
        cell_table.FirstSigBin(idx(cIdx),:) = first_sig_bin;
        
        % get peak response       
        psth_resp_laser_this = squeeze(shuf_test.psth_laser(cIdx,:,t_psth>0 & t_psth<=pulseDurMs));
        psth_resp_laser_this = psth_resp_laser_this(pulse_reorder,:);
        peak_resp_bin = nan(nCond * nPulse,1);
        for condIdx = 1:numel(peak_resp_bin)
            [~,peak_resp_bin(condIdx)] = max(psth_resp_laser_this(condIdx,:));
        end
        cell_table.PeakResponseBin(idx(cIdx),:) = peak_resp_bin;
        

        responsive_pulse = find(shuf_test.sig_laser_response(cIdx,pulse_reorder));

        trigger_collision = 1000/opt.bin_size_ms * laser_ts_all';
        grp = reshape(repmat((TrialTypes-1)*nPulse,nPulse,1)+repmat((1:nPulse)',1,numel(TrialTypes)),numel(trigger_collision),1);

        keep = ismember(grp,responsive_pulse);
        trigger_collision = trigger_collision(keep);
        grp = grp(keep);

        trigger_collision = trigger_collision + first_sig_bin(grp)/shuf_test.opt.bin_size_ms;
        grp = ceil(grp/10);

        [~,psth] = plot_timecourse('timestamp',spiket,trigger_collision,trigger_collision-10,trigger_collision+10,grp,...
            'win_len',1,'resample_bin',1,'plot_type','none');
        
        if ~isempty(psth)
            window = psth.x*opt.bin_size_ms*1000 >= opt.collision_test_window(1) & ...
                psth.x*opt.bin_size_ms*1000 <= opt.collision_test_window(2);
            trials_with_spikes = sum(psth.rate_rsp(:,window),2)>0;

            grp2 = grp;
            grp2(trials_with_spikes) = 1;
            grp2(~trials_with_spikes) = 0;

            [~,psth] = plot_timecourse('timestamp',spiket,trigger_collision,trigger_collision-10,trigger_collision+10,grp2,...
                'win_len',1,'resample_bin',1,'plot_type','none');

            if size(psth.mean,1)==2
                collision_test(cIdx) = psth.mean(2,psth.x==0) == 0;
            else
                collision_test(cIdx) = 1;
            end
        else
            collision_test(cIdx) = 0;
        end
        
        % IS THE CELL OPTOTAGGED?
        % 1) QC
        pass_qc = metrics.isi_viol(cIdx) < opt.isi_viol_thresh & ...
            metrics.amplitude_cutoff(cIdx) < opt.amplitude_cutoff_thresh & ...
            metrics.presence_ratio(cIdx) > opt.presence_ratio_thresh;

        % 2) laser responsive and pass firing rate cutoff (during laser pulse)
        laser_responsive = any(laser_response(cIdx,:) & pass_fr_cutoff(cIdx,:));

        % 3) waveform correlation
        pass_wv_corr = wv_corr.waveform_corr(cIdx) > opt.waveform_corr_thresh;

        % 4) collision test; optotagged units pass all tests:
        optotagged(cIdx) = pass_qc & laser_responsive & pass_wv_corr & collision_test(cIdx);
        
        % Make figs
        if opt.make_figs && optotagged(cIdx)
            hfig = figure('Position',[200 200 1200 400]);
            hfig.Name = sprintf('%s_c%d',opt.session,good_cells(cIdx));

            % laser psth
            subplot(1,3,1);
            ax = plot_timecourse('timestamp',spiket,trigger,trigger-opt.ms_before/opt.bin_size_ms,trigger+opt.ms_after/opt.bin_size_ms,...
                TrialTypes','win_len',5,'resample_bin',5,'plot_type','both');

            % make figure pretty
            ax(1).Title.String = sprintf('LASER_RESP=[%d %d %d %d]',laser_response(cIdx,:) & pass_fr_cutoff(cIdx,:));
            ax(1).Title.Interpreter = 'none';
            ax(2).Legend.String = laser_target(laser_reorder);
            xtick = -opt.ms_before:100:opt.ms_after;
            ax(2).XTick = xtick/1000/opt.bin_size_ms;
            xticklabel = cell(numel(xtick));
            for tickNum = 1:numel(xtick)
                xticklabel{tickNum} = num2str(xtick(tickNum));
            end
            ax(2).XTickLabel = xticklabel;
            ax(2).XLabel.String = 'ms from laser onset';

            % waveform corr
            subplot(1,3,2); hold on;
            plot(wv_corr.non_laser_evoked_waveforms(cIdx,:));
            plot(wv_corr.laser_evoked_waveforms(cIdx,:));
            title(sprintf('%s_c%d, ISI_VIOL=%0.2f, PRES_RATIO=%0.2f, AMP_CUTOFF=%0.2f, OPTOTAG=%d\nWV_CORR=%0.3f',...
                opt.session,good_cells(cIdx),metrics.isi_viol(cIdx),metrics.presence_ratio(cIdx),metrics.amplitude_cutoff(cIdx),...
                optotagged(cIdx),wv_corr.waveform_corr(cIdx)),'Interpreter','none');
            legend({'non-laser','laser','all'},'Location','southeast');
            xlabel('samples');
            ylabel('uV');
            xlim([1 size(wv_corr.laser_evoked_waveforms,2)]);

            % collision test
            subplot(1,3,3);
            if ~isempty(trigger_collision)
                [ax,psth] = plot_timecourse('timestamp',spiket,trigger_collision,trigger_collision-10,trigger_collision+10,grp2,...
                    'win_len',1,'resample_bin',1,'plot_type','both');
                ax(1).Title.String = sprintf('COLLISION_TEST=%d',collision_test(cIdx));
                ax(1).Title.Interpreter = 'none';
                if size(psth.mean,1)==2
                    ax(2).Legend.String = {'no collision','collision'};
                end
            else
                title(sprintf('COLLISION_TEST=%d',collision_test(cIdx)),'Interpreter','none');
                xticks([]);
                yticks([]);
            end
            
            saveas(hfig,fullfile(paths.figs,hfig.Name),'png');
            close(hfig);
            
        end

    end
    cell_table.CollisionTest(idx) = collision_test;
    cell_table.Optotagged(idx) = optotagged;
    
end
cell_table = cell_table(1:ctr-1,:);
toc

%% test criteria

% 1) QC
pass_qc = cell_table.isi_viol < opt.isi_viol_thresh & ...
    cell_table.amplitude_cutoff < opt.amplitude_cutoff_thresh & ...
    cell_table.presence_ratio > opt.presence_ratio_thresh;

% 2) laser responsive and pass firing rate cutoff (during laser pulse)
laser_responsive = any(cell_table.LaserResponse & cell_table.PassFRCutoff,2);

% 3) waveform correlation
pass_wv_corr = cell_table.waveform_corr > opt.waveform_corr_thresh;

% 4) collision test
collision_test = cell_table.CollisionTest;

% optotagged units pass all tests:
optotagged = pass_qc & laser_responsive & pass_wv_corr & collision_test;

cell_table_all = cell_table;
cell_table_qc = cell_table(pass_qc,:);
cell_table_opto = cell_table(optotagged,:);
