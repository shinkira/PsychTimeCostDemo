classdef publish_plot < handle
    properties
        h_ax
        h_fig
        figure_label
        active_axis
        data
        legend
        text
        savename
        savename_fig
        savedir
        invert_colors = 0;
    end
    
    methods
        function obj = publish_plot(n_filas,n_col,varargin)
            hfig = [];
            for i=1:length(varargin)
                if isequal(varargin{i},'hfig')
                    hfig = varargin{i+1};
                end
            end
            
            if ~isempty(hfig)
                obj.h_fig = figure(hfig);
            else
                obj.h_fig = figure;
            end
            
            cont = 0;
            for i = 1:n_filas
                for j = 1:n_col
                    cont = cont + 1;
                    obj.h_ax(cont) = subplot(n_filas,n_col,cont);
                end
            end
            obj.data.n_filas = n_filas;
            obj.data.n_col = n_col;
            
            set(gcf,'color','w')
        end
        
        
        
        function load_from_fig(obj,figure_filename)
            h_fig_old = obj.h_fig;
            h_fig   = open(figure_filename);
            children = get(h_fig,'children');
            h_ax    = children(end:-1:1);
            obj.h_fig = h_fig;
            obj.h_ax  = h_ax;
            set(obj.h_fig,'color','w')
            close(h_fig_old)
            
        end
        
        function set_fig_size(obj,varargin)
            
            xSize = 21.6; ySize = 27.9;%tama�o de la figure
            for i = 1:length(varargin)
                if isequal(varargin{i},'xSize')
                    xSize = varargin{i+1};
                elseif isequal(varargin{i},'ySize')
                    ySize = varargin{i+1};
                end
            end
            
            % Tama�o y posici�n de la figura
            set(obj.h_fig,'PaperUnits','centimeters')
            
            xLeft = (21.6-xSize)/2; yTop = (27.9-ySize)/2;
            set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
            
            %Para lo que se ve en pantalla:
            set(gcf,'Position',[200    -100   700*xSize/21.6  888*ySize/27.9])
            
        end
        
        function new_axes(obj,varargin)
            position = [];
            for i = 1:length(varargin)
                if isequal(varargin{i},'position') || isequal(varargin{i},'Position')
                    position = varargin{i+1};
                end
            end
            if isempty(position)
                obj.h_ax(end+1) = axes();
            else
                obj.h_ax(end+1) = axes('position',position);
            end
        end
        
        function resize_vertical(obj,ax_num,vert_ratio,vert_separation)
            % obj.resize_vertical(obj,[1,2],[1,3])
            
            % vertical separation
            if nargin<4
                vert_separation = 0.13;
            end
            
            % get positions
            pos = [];
            for i=1:length(ax_num)
                pos = [pos; get(obj.h_ax(ax_num(i)),'Position')];
            end
            % resort
            [~,ind] = sort(pos(:,2));
            pos = pos(ind,:);
            ax_num = ax_num(ind);
            
            % new vert sizes
            V = pos(:,4);
            newV = V(:).*vert_ratio(:)/sum(vert_ratio);
            newV = newV * sum(V)/sum(newV);
            
            
            
            % new positions
            for i=1:length(ax_num)
                posit(i,:) = [pos(i,1),pos(i,2),pos(i,3),newV(i)];
                if i>1
                    posit(i,2) = posit(i-1,2)+posit(i-1,4)+vert_separation;
                end
            end
            % go
            for i=1:length(ax_num)
                set(obj.h_ax(ax_num(i)),'Position',posit(i,:));
            end
        end
        
        
        function align(obj)
            % align top, bottom, left, right
            % to do
        end
        function distribute(obj)
            % to do
        end
        
        function draggable(obj,h)
            
            right = 29;
            left  = 28;
            up    = 30;
            down  = 31;
            DELTA = [0.05 0.01];
            delta = DELTA(1);
            
            key_d = 100; % letra d: duplica axes
            key_n = 110;
            key_w = 119;
            
            var      = 'position';
            %             var      = 'size';
            char  = getkey(1);
            position = get(h,'position');
            figure(obj.h_fig)
            h_borde = [];
            while not(char==27)
                if char==key_n % next axes
                    obj.next();
                    h = gca;
                    position = get(h,'position');
                end
                if char==key_d
                    % duplicates axes
                    obj.new_axes('position',position);
                    h = obj.h_ax(end);
                    position = get(h,'position');
                end
                if char==key_w
                    ind = find(DELTA == delta);
                    if ind==length(DELTA)
                        delta = DELTA(1);
                    else
                        delta = DELTA(ind+1);
                    end
                end
                
                switch char
                    case 115 %"s" for size
                        var = 'size';
                        char  = getkey(1);
                    case 112 % "p" for position
                        var = 'position';
                        char  = getkey(1);
                end
                if isequal(var,'position')
                    switch char
                        case right
                            position(1) = position(1) + delta;
                        case up
                            position(2) = position(2) + delta;
                        case left
                            position(1) = position(1) - delta;
                        case down
                            position(2) = position(2) - delta;
                    end
                elseif isequal(var,'size')
                    switch char
                        case right
                            position(3) = position(3) + delta;
                        case up
                            position(4) = position(4) + delta;
                        case left
                            position(3) = position(3) - delta;
                        case down
                            position(4) = position(4) - delta;
                    end
                end
                
                set(h,'position',position)
                delete(h_borde);
                h_borde = borde_en_axis('ax',h);
                position = get(h,'position');
                
                char  = getkey(1);
                figure(obj.h_fig)
            end
            delete(h_borde);
            
        end
        
        function draggable_print(obj)
            for i=1:length(obj.h_ax)
                pos = get(obj.h_ax(i),'Position');
                str= ['set(p.h_ax(',num2str(i),'),''Position'',[',num2str(pos),'])'];
                disp(str)
            end
        end
        
        function copy_from_ax(obj,h_from,h_to)
            if nargin<3 || isempty(h_to)
                h_to = obj.h_ax(obj.active_axis);
            end
            copyobj(allchild(h_from),h_to);
            %             copyobj(h_from,h_to);
            props = {'xlim','ylim','xlabel','ylabel','yaxislocation'};
            for i=1:length(props)
                if isempty(strfind(version,'R2014b')) && (isequal(props{i},'xlabel') || isequal(props{i},'ylabel'))
                    gg = get(get(h_from,props{i}),'String');
                    set(get(h_from,props{i}),'String',gg);
                else
                    prop = get(h_from,props{i});
                    set(h_to,props{i},prop);
                end
                
            end
            %             xli = get(h_from,'xlim');
            %             yli = get(h_from,'ylim');
            %             set(h_to,'xlim',xli);
            %             set(h_to,'ylim',yli);
            
        end
        
        function saveas(obj,varargin)
            % save as .fig
            
            filename_fig = 'f.fig';
            if length(varargin)==1
                filename_fig = varargin{1};
            else
                for i = 1:length(varargin)
                    if isequal(varargin{i},'filename')
                        filename_fig = varargin{i+1};
                    end
                end
            end
            
            obj.data.filename_fig = filename_fig;
            saveas(obj.h_fig,filename_fig);
            obj.savename_fig = filename_fig;
            
            
        end
        
        
        function plot(obj,i,varargin)
            set(gcf,'CurrentAxes',obj.h_ax(i))
            plot(varargin{:},'LineWidth',0.5)
            set(obj.h_ax(i),'FontSize',6.8,'LineWidth',0.335)
            set(gca,'box','off')
            %obj.data.data(i) = varargin;
            %obj.data.fun(i) = 'plot';
        end
        
        function set_active(obj,varargin)
            obj.current_ax(obj,varargin{2:end});
        end
        
        function current_ax(obj,varargin)
            if length(varargin)==1
                i = varargin{1};
                set(obj.h_fig,'CurrentAxes',obj.h_ax(i))
                %                 obj.active_axis = obj.h_ax(i);
                obj.active_axis = i;
            elseif length(varargin)==2
                i = varargin{1}; j = varargin{2};
                z = obj.data.n_col*(i-1)+j;
                set(obj.h_fig,'CurrentAxes',obj.h_ax(z))
                %                 obj.active_axis = obj.h_ax(z);
                obj.active_axis = z;
            else
                disp('error en input')
            end
        end
        
        function delete_ax(obj,h)
            if mod(h,1)==0
                ind = h;
                h   = obj.h_ax(ind);
            else
                ind = find(obj.h_ax==h);
            end
            delete(h); % deletes the axes
            obj.h_ax(ind) = [];
            
        end
        
        function next(obj)
            if isempty(obj.active_axis)
                set(obj.h_fig,'CurrentAxes',obj.h_ax(1))
                obj.active_axis = 1;
            else
                obj.active_axis = obj.active_axis+1;
                if obj.active_axis>length(obj.h_ax)
                    obj.active_axis = 1;
                end
                set(obj.h_fig,'CurrentAxes',obj.h_ax(obj.active_axis))
            end
        end
        
        function save(obj,varargin)
            %p.save('path',pwd,'ancho',7,'filename','blabla','renderer','painters','resolution',150)
            filename = 'aaa';
            dire = pwd;
            renderer = 'painters';
            %             resolution = num2str(150);
            %             ancho = 21.6;
            set_ancho = false;
            for i = 1:length(varargin)
                switch varargin{i}
                    case 'filename'
                        filename = varargin{i+1};
                    case 'path'
                        dire = varargin{i+1};
                    case 'ancho'
                        ancho  = varargin{i+1};
                        set_ancho = false;
                    case 'renderer'
                        renderer = varargin{i+1};
                        %                     case 'resolution'
                        %                         resolution = num2str(varargin{i+1});
                end
            end
            if set_ancho
                set_figure_size(ancho);
            end
            
            obj.savename = filename;
            obj.savedir = dire;
            old_dir = pwd;
            cd(dire)
            figName = [filename,'.eps'];
            % eval(['print -depsc2 -tiff -',renderer,' -r',resolution,' ',figName])
            eval(['print -depsc2 -tiff -',renderer,' ',figName])
            cd(old_dir)
        end
        
        function format(obj,varargin)
            FontSize = 6.8;
            LineWidthAxes = 0.335;
            LineWidthPlot = [];
            for i=1:length(varargin)
                if isequal(varargin{i},'presentation')
                    FontSize = 20;
                    LineWidthAxes = 1;
                    LineWidthPlot = 2.5;
                elseif isequal(varargin{i},'invert_colors') && obj.invert_colors==0
                    obj.invert_colors = 1;
                elseif isequal(varargin{i},'FontSize')
                    FontSize = varargin{i+1};
                elseif isequal(varargin{i},'LineWidthPlot')
                    LineWidthPlot = varargin{i+1};
                elseif isequal(varargin{i},'LineWidthAxes')
                    LineWidthAxes = varargin{i+1};
                end
            end
            
            if not(isempty(LineWidthPlot))
                for j=1:length(obj.h_ax)
                    set(get(obj.h_ax(j),'Children'),'LineWidth',LineWidthPlot)
                end
            end
            
            all_text = findall(obj.h_fig,'Type','text');
            set(all_text,'FontSize',FontSize)
            
            for i=1:length(obj.h_ax)
                set(obj.h_ax(i),'FontSize',FontSize,'LineWidth',LineWidthAxes,'box','off')
            end
            
            
            if obj.invert_colors==1
                a = findall(obj.h_fig);
                w = findobj(a,'Color','w');
                b = findobj(a,'Color','k');
                set(w,'Color','k');
                set(b,'Color','w');
                
                for j=1:length(obj.h_ax)
                    set(obj.h_ax(j),'Ycolor','w')
                    set(obj.h_ax(j),'Xcolor','w')
                    set(obj.h_ax(j),'Color','none') % Para fondo
                    %transparente
                end
                
                set(obj.h_fig,'Color','none') % Para fondo transparente
                set(obj.h_fig,'InvertHardcopy','off')
            end
            
        end
        
        function number_the_plots(obj,action)
            if isequal(action,'show')
                obj.text.number_plot = [];
                for i = 1:length(obj.h_ax)
                    obj.current_ax(i);
                    xli = xlim;
                    yli = ylim;
                    obj.text.number_plot(i) = text(sum(xli)/2,sum(yli)/2,num2str(i));
                end
            elseif isequal(action,'hide')
                delete(obj.text.number_plot);
            end
            
            
            
        end
        
        function shadow_plot(obj,i,varargin)
            set(gcf,'CurrentAxes',obj.h_ax(i))
            [errorPatch,dataLine] = niceBars(varargin{:});
            set(obj.h_ax(i),'FontSize',6.8,'LineWidth',0.335)
            set(gca,'box','off')
            obj.data(i).data = varargin;
            obj.data(i).fun = 'niceBars';
        end
        
        function rainbow_plot(obj,i,x,y,varargin)
            set(gcf,'CurrentAxes',obj.h_ax(i))
            rplot(x,y,varargin{:});
            set(obj.h_ax(i),'FontSize',6.8,'LineWidth',0.335)
            set(gca,'box','off')
            obj.data(i).data = {x,y,varargin{:}};
            obj.data(i).fun = 'rainbow_plot';
        end
        
        function legend_save(obj)
            fid = fopen([obj.savedir,'/',['legend_',obj.savename,'.txt']],'w');
            fprintf(fid,'%s\n','\begin{figure}[t]');
            fprintf(fid,'%s\n','\begin{center}');
            fprintf(fid,'%s\n',['\includegraphics{',obj.savename,'.eps}']);
            fprintf(fid,'%s\n',['\caption{',obj.legend,'}']);
            fprintf(fid,'%s\n',['\label{',obj.savename,'}']);
            fprintf(fid,'%s\n','\end{center}');
            fprintf(fid,'%s\n','\end{figure}');
            fclose(fid);
        end
        
        function displace_ax(obj,h_ax,delta_pos)
            
            if all(mod(h_ax,1)==0)
                h_ax = obj.h_ax(h_ax);
            end
            
            for i=1:length(h_ax)
                pos = get(h_ax(i),'Position');
                newpos = pos + delta_pos;
                set(h_ax(i),'Position',newpos);
            end
            
        end
        
        function hax = handle(obj,i,j)
            ind = obj.data.n_col*(i-1)+j;
            hax = obj.h_ax(ind);
            
        end
        
        
        function plot_for_presentation
        end
        
        function save_for_presentation
        end
    end
end
