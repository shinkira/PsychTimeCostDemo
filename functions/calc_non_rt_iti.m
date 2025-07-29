function [non_rt_iti_val, non_rt_iti_inv, non_rt_iti_brk] = calc_non_rt_iti(data)

    % Copied from psychIntervalCalc
    % non_rt_val & non_rt_val include inter-trial intervals (~2sec feedback)

    % data = [result', response', coh', rt', start_t', end_t', pts', max_dot_dur']

    % R(strfind(str_result,'WRONG'))    = 0;
    % R(strfind(str_result,'CORRECT'))  = 1;
    % R(strfind(str_result,'NOFIX'))    = 2;
    % R(strfind(str_result,'FIXBREAK')) = 3;
    % R(strfind(str_result,'NOCHOICE')) = 4;
    % R(strfind(str_result,'TIMEOUT'))  = 5;

    % Calculated in "sec"

%     val_ind      = data(:,1)==0 | data(:,1)==1;
%     time_lim_ind = data(:,1)==5;
%     break_ind    = data(:,1)==2 | data(:,1)==3 | data(:,1)==4;
    
    val_ind      = ismember(data(:,1)',[0,1]);
    time_lim_ind = ismember(data(:,1)',[5]);
    break_ind    = ismember(data(:,1)',[2,3,4]);

    trial_length = data(:,6) - data(:,5);
    rt = data(:,4)/1000;
    max_dot_dur = data(:,8);
    iti = data(2:end,5) - data(1:(end-1),6);
    iti = iti(iti>1);
    iti_mean = mean(iti);

    non_rt_val = trial_length(val_ind) - rt(val_ind);
    non_rt_inv = trial_length(time_lim_ind) - max_dot_dur(time_lim_ind);
    non_rt_brk = trial_length(break_ind) - rt(break_ind);
    
    non_rt_iti_val = mean(non_rt_val,'omitnan') + iti_mean;
    non_rt_iti_inv = mean(non_rt_inv,'omitnan') + iti_mean;
    non_rt_iti_brk = mean(non_rt_brk,'omitnan') + iti_mean;
    
end 