function fit_nDDM(debug)
    
    if ~exist('debug','var') || isempty(debug)
        debug = false;
    end

    p = mfilename('fullpath');
    [p,~,~] = fileparts(p);
    cd(p);

    % non-parametric DDM (nDDM) has three fitted parameters in theta.
    % theta(1): kappa
    % theta(2): mean of non-decision time
    % theta(3): s.d. of non-decision time

    % tl: lowest value for theata for fit
    % th: highest value for theta for fit

    % Initiall guesses for fit are drawn from a uniform distribution between 
    % the lower bound (tl) and upper bound (th).
    % During the fit, tl and th also provides the lower and upper bounds.
    % If they are the same value, the parameter is fixed at that value during 
    % the fit. See fminsearchbnd.m for details.

    tl = [  1, 0.001, 0.001, 0, 0];
    th = [100,     1,     1, 0, 0];
    % tl = [ 5, 0.1 ,0.001, 0, 0];
    % th = [40, 0.4 ,0.08,  0, 0];
    tg = [];

    fn_wrapper = @wrapper_dtb_empiricalbound_rt;
    fit_diffusion(fn_wrapper,tl,th,tg,debug);

end
