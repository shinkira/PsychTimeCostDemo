function P = dtb_fp_cn_vec(drift,t,Bup,Blo,y,y0,notabs_flag,varargin)
% crank nicholson, 1d advection diffusion

% https://georg.io/2013/12/03/Crank_Nicolson
% "The Crank-Nicolson method implemented from scratch in Python"

% Expired link: http://georg.io/2013/12/Crank_Nicolson_Convection_Diffusion/ 

% First arg in varargin is demo_flag

verbose = 0;

if isempty(varargin)
    demo_flag = 0;
else
    demo_flag = varargin{1};
end

win = 4;
if demo_flag
    t = movmean(t,win+1,'Endpoints','discard');
    t = t(1:win:end);
    Bup = movmean(Bup,win+1,'Endpoints','discard');
    Bup = Bup(1:win:end);
    Blo = movmean(Blo,win+1,'Endpoints','discard');
    Blo = Blo(1:win:end);
    
    win_y = 2;
    y  =  y(1:win_y:end);
    y0 = y0(1:win_y:end);
end

dt = t(2) - t(1);
dy = y(2) - y(1);
ny = length(y);
nt = length(t);


% Expand any flat bounds
if numel(Bup)==1
    Bup = repmat(Bup,nt,1);
end
if numel(Blo)==1
    Blo = repmat(Blo,nt,1);
end

P = struct('drift',drift,'t',t,'Bup',Bup,'Blo',Blo,'y',y,'y0',y0,...
    'notabs_flag',notabs_flag);

nd = length(drift);

% Preallocate
P.up.pdf_t = zeros(nd,nt);
P.lo.pdf_t = zeros(nd,nt);
if notabs_flag
    P.notabs.pdf = zeros(nd,ny,nt);
end
p_threshold = 1.0E-5; % Threshold for proportion un-terminated to stop simulation

% params advective diffusion equation

[A,B] = crank_nicholson_sparsematrix(drift,nd,ny,dy,dt);

yr = repmat(y(:),nd,1);

u = repmat(y0(:),nd,1);

for k = 2:nt
    
    d = B*u;
    u = A\d;
    
    ur = reshape(u,ny,nd); % maybe too slow?
    
    % Select density that has crossed bounds
    P.up.pdf_t(:,k) = sum(ur(y>=Bup(k),:),1);
    P.lo.pdf_t(:,k) = sum(ur(y<=Blo(k),:),1);

    % Keep only density within bounds
    outofbounds = yr<=Blo(k) | yr>=Bup(k);
    u(outofbounds) = 0;

    % Save if requested
    if notabs_flag
        P.notabs.pdf(:,:,k) = ur';
        ur_cut = ur;
        ur_cut(y>=Bup(k) | y<=Blo(k),:) = 0;
    end
    
    if sum(sum(ur,1)<p_threshold)==nd
        break;
    end
    
    if verbose && mod(k,1e2)==0
        fprintf('%d / %d\n',k,nt)
    end
end


if notabs_flag
    P.notabs.pos_t = sum(P.notabs.pdf(:,y'>=0,:),2);
    P.notabs.neg_t = sum(P.notabs.pdf(:,y'< 0,:),2);
end

P.up.p = sum(P.up.pdf_t,2);
P.lo.p = sum(P.lo.pdf_t,2);

t = t(:);

P.up.mean_t = transpose(t'*P.up.pdf_t')./P.up.p;
P.lo.mean_t = transpose(t'*P.lo.pdf_t')./P.lo.p;

P.up.cdf_t = cumsum(P.up.pdf_t,2);
P.lo.cdf_t = cumsum(P.lo.pdf_t,2);
