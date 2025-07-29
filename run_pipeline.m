clear
debug = 0;

%% Fit non-parametric Drift Diffusion Model (nDDM) to behavioral data of subjects
fit_nDDM(debug);

%% Compile MEX code for Dynamic Programming
% You need to install MATLAB Coder: https://www.mathworks.com/help/coder/
% In command menu, click Add-Ons > Get Add-Ons
compile_d2b_dp;

%% Using the fitted parameters in nDDM, find the optimal bound shape with Dyanamic Programming
run_optim_bound_dp(debug);

%% Retrieve and save the bounds derived by nDDM and Dynamic Programming 
plotBound(1,debug);

%% Plot the bounds derived by nDDM and Dynamic Programming 
plotBound(2,debug);