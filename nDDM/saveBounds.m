function saveBounds

    vSubject = {'A','B','C','N'};
    vStage = {'Pre','TimeLimit'};
    IND = combvec(1:length(vSubject),1:length(vStage));

    subjects = vSubject(IND(:,1));
    stages   = vStage(IND(:,2));
    
    [p,~,~] = fileparts(mfilename('fullpath'));
    data_dir = fullfile(p,'data');
    
    for si = 1:4
        for ti = 1:2
            file_name = sprintf('%s_%s.mat',vSubject{si},vStage{ti});
            load(fullfile(data_dir,file_name),'P');
            B.(vSubject{si}).(vStage{ti}) = P.Bup;
            Theta.(vSubject{si}).(vStage{ti}) = P.optim.theta;
            % B.(subjectsP1{g.group(i)}).Pre;
        end
    end
    save(fullfile(data_dir,'EmpiricalBound.mat'),'B','Theta');
end