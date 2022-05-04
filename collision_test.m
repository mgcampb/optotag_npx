function [test_res, pval] = collision_test(spiket, trigger, rsp_win, coll_win, varargin)
% Performs collision test 
% MGC 4/14/2022
%
%   INPUT:
%       spiket      :   Spike times in ms
%       trigger     :   Trigger times in ms
%       resp_win    :   Response window in ms [winstart winend]
%       coll_win    :   Collision window in ms [winstart winend]
%
%   OUTPUT:
%       test_res    :   1 if passes, 0 if fails
%       pval        :   p-value of rank sum test between spike counts in
%                       response window for trials with spikes in collision
%                       window versus not

    make_plot = false;
    process_varargin(varargin);
    
    has_coll_spk = zeros(numel(trigger),1);
    for i = 1:numel(trigger)
        if sum(spiket>trigger(i)+coll_win(1) & spiket<trigger(i)+coll_win(2)) > 0 
            has_coll_spk(i) = 1;
        end
    end

    num_rsp_spk = nan(numel(trigger),1);
    for i = 1:numel(trigger)
        num_rsp_spk(i) = sum(spiket>trigger(i)+rsp_win(1) & spiket<trigger(i)+rsp_win(2));
    end

    if sum(has_coll_spk) == 0
        test_res = 1;
    else % possibly make this less strict
        if sum(num_rsp_spk(has_coll_spk==1)) == 0
            test_res = 1;
        else
            test_res = 0;
        end
    end

    pval = ranksum(num_rsp_spk(has_coll_spk==1),num_rsp_spk(has_coll_spk==0));

    if make_plot
        figure;
        plot_timecourse('timestamp',spiket,trigger+rsp_win(1),-10,10,has_coll_spk,'win_len',1,'resample_bin',1);
    end

end