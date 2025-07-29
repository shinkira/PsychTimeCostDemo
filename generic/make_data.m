function dat = make_data(subject,stage,flag,pars)

if nargin<3 || strcmp(stage,'Pre') || strcmp(stage,'Post')
    flag = [];
end

load_dir_name = fullfile(pars.root,'SubjectData');
save_dir_name = fullfile(pars.root,'MiniData');
makedir(save_dir_name);

% dir_name = '../../generic/data';
if isempty(flag)
    file_name = sprintf('%s_%s.mat',subject,stage);
else
    file_name = sprintf('%s_%s_%s.mat',subject,stage,flag);
end

if ~exist(fullfile(save_dir_name,file_name),'file')

    load(fullfile(load_dir_name,'psychInfoNoDots.mat'))

    %% data
    % 1: result (1:correct, 0 error)
    % 2: choice direction (1:right, 2:left)
    % 3: coherence (pos: right, neg:left)
    % 4: RT (ms)
    % 5: time of trial start
    % 6: time of trial end
    % 7: point (start from 50)
    % 8: deadline (s)
    
    % R(strfind(str_result,'WRONG'))    = 0;
    % R(strfind(str_result,'CORRECT'))  = 1;
    % R(strfind(str_result,'NOFIX'))    = 2;
    % R(strfind(str_result,'FIXBREAK')) = 3;
    % R(strfind(str_result,'NOCHOICE')) = 4;
    % R(strfind(str_result,'TIMEOUT'))  = 5;
    
    if isempty(flag)
        
        data = (psychInfo.(subject).(stage).N.selected_data);
        choice = data(:,2)==1;
        coh = data(:,3)/1000;
        c = data(:,1);
        rt = data(:,4)/1000;

        inds = ismember(c,[0,1]) & ~isnan(rt);
        dat.choice = choice(inds);
        dat.coh    = coh(inds);
        dat.rt     = rt(inds);
        dat.c      = c(inds);
        
    elseif strcmp(flag,'all')
        
        data = (psychInfo.(subject).(stage).N.data); % include time limit trials
        choice = data(:,2)==1;
        coh = data(:,3)/1000;
        c = data(:,1);
        rt = data(:,4)/1000;
        pt = data(:,7);
        dl = data(:,8);

        dat.choice = choice;
        dat.coh    = coh;
        dat.rt     = rt;
        dat.c      = c;
        dat.pt     = pt;
        dat.dl     = dl;
        
    end

    save(fullfile(save_dir_name,file_name),'dat');

else
    
    load(fullfile(save_dir_name,file_name),'dat');    
    
end

