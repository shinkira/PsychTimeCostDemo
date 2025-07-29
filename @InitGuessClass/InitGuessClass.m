classdef InitGuessClass < handle
    properties
        tl
        th
        tg
        TG
        tguess_flag
        group = [];
        N
        temp_filename
        gather_filename
        pre_filename = 'optim'
    end
    
    methods
        function obj = InitGuessClass(tl,th,tg)
            if nargin<3
                tg = [];
            end
            obj.tl = tl;
            obj.th = th;
            obj.tg = tg;
            obj.N = 1; %added

        end
        
        function generate_tguess_current(obj,uni_group)
            if nargin<2
                uni_group = 1;
            end
            if isscalar(uni_group)
                uni_group = 1:uni_group;
            end
            
            obj.TG = [];
            obj.group = [];
            for i=1:length(uni_group)
                obj.TG(i,:) = obj.tg;
                obj.group(i) = uni_group(i);
            end
            obj.N = length(obj.group);
            obj.generate_randomfilenames();
        end
            
        function generate_tguess_random(obj,nguesses,seed,uni_group)
            if nargin<4
                uni_group = 1;
            end
            if isscalar(uni_group)
                uni_group = 1:uni_group;
            end
            
            %set seed
            rstream = RandStream('mt19937ar','Seed',seed);
            
            obj.TG = [];
            obj.group = [];
            cont = 0;
            for i=1:length(uni_group)
                for j=1:nguesses
                    cont = cont + 1;
                    obj.TG(cont,:) = obj.sample_tguess(obj.tl,obj.th,rstream);
                    % obj.TG(cont,:) = diffusion_fun.sample_tguess(obj.tl,obj.th,rstream);
                    obj.group(cont) = uni_group(i);
                end
            end
            
            obj.N = length(obj.group);
            obj.generate_randomfilenames();
        end
        
        
        function keep_best(obj,fn_handle,Nkeep)
            n = obj.N;
            for i=1:n
                val(i) = feval(fn_handle,obj.TG(i,:));
            end
            [~,ind_keep] = sort(val);
            ind_keep = ind_keep(1:Nkeep);
            
            % replace with the ones to keep
            obj.TG = obj.TG(ind_keep,:);
            obj.group = obj.group(ind_keep);
            obj.N = Nkeep;
            obj.temp_filename = obj.temp_filename(ind_keep);
            
        end
        
        
        function generate_randomfilenames(obj)
            % generates random temp filenames
            for i=1:obj.N
                temp_filename = [obj.pre_filename,'_group',num2str(obj.group(i)),'_r',num2str(ceil(rand*intmax))];
                obj.temp_filename{i} = temp_filename;
            end
        end
        
        function get_optim_filename(obj,n)
            for i=1:n
                savefilename = [obj.pre_filename,'_group',num2str(i)];
                obj.gather_filename{i} = savefilename;
            end
        end
        
        function gather_files(obj)
            %get the min fval for each group
            
            for i=1:obj.N
                f = obj.temp_filename{i};
                aux = load(f);
                fnames = fieldnames(aux);
                for j=1:length(fnames)
                    all(i).(fnames{j}) = aux.(fnames{j});
                    all(i).group = obj.group(i);
                end
            end
            
            % best per group
            fval = [all.fval];
            u = unique(obj.group);
            for i=1:length(u)
                I = find(obj.group==u(i));
                [~,ind] = min(fval(I));
                best_id = I(ind);
                todos = all(I);
                mejor = all(best_id);
                savefilename = [obj.pre_filename,'_group',num2str(i)];
                obj.gather_filename{i} = savefilename;
                save(savefilename,'todos','mejor');
            end
            
        end
        
        function delete_temp_files(obj)
            for i=1:obj.N
                filename = obj.temp_filename{i};
                if isunix
                    eval(['!rm ',filename,'.mat']);
                else
                    % for MATALB
                    delete([filename,'.mat']);
                end
            end
        end
        
        function delete_gather_files(obj)
            for i=1:length(unique(obj.group))
                filename = obj.gather_filename{i};
                eval(['!rm ',filename,'.mat']);
            end
        end
        
        function generate_tguess_grid(obj)
            % to do
            
        end
    end
end