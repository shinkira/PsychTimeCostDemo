function [bestQ,Q,V,Vdelta,rho,TlistenVal] = d2b_dp(MAXDV,DVDELTA,maxT,dt,coherences,kappa,gamma,rewardPrm,maxError,p_dl,ITI,policy,TlistenVal2D)

    % INPUT ARGUMENTS:
    % MAXDV:        The range of decision variable [-MAXDV, MAXDV] used to compute
    %               the optimal bound shape. Set it higher than what the actual bound
    %               height ought to be.
    % 
    % DVDELTA:      Grid size for the deccision variable
    % 
    % maxT:         The range of trial length [0, maxT] used to compute the optimal
    %               bound shape. Set it longer than the maximum RT.
    % 
    % dt:           Grid size for the elapsed time in a trial
    % 
    % coherences:   Signed coherences of the random dot stimulus (in fraction, not percent)
    % 
    % kappa:        Scaling factor for the drift term
    % 
    % gamma:        Temporal discounting factor. Set to zero if it is not discounted.
    % 
    % rewardPrm:    Reward associated with the correct choice, error choice,
    %               fixation, and time-up
    % 
    % maxError:     Maximum error allowed for the estimation of the reward rate
    % 
    % p_dl:         Probability of the deadline occuring at a given moment when it
    %               has not occured yet (i.e., hazard rate of the provisional deadline)
    % 
    % ITI:          Mean inter-trial interval defined as the time between the offset 
    %               and onset of the random-dot stimulus (including the
    %               time to acquire fixation and the stimulus latencies)
    %               ITI(1) is for valid trials and ITI(2) ismin for canceled
    %               trials (ITI_inv)
    % 
    % policy:       A 2D matrix (like bestQ) specifying the action to take in each sate. 
    %               if provieded, the function will compute the reward rate for this policy
    % 
    % TlistenVal2D: Transition matrix 
    %
    %
    % OUTPUT ARGUMENTS:
    % bestQ:        A 2D matrix stroing the best action in a given state
    % 
    % Q:            A 3D matrix storing the action values for three actions in a given sate
    % 
    % V:            A 2D matrix of the state value for a given state
    % 
    % rho:          Reward rate (amount of reward obtained per second)
    % 
    % TlistenVal:   A 3D transition matrix, where TlistenVal(i,j,k) is the 
    %               transition probability of the decision variable (DV) from DV(i) to DV(k) 
    %               at time step t(j).
    % 
    % if TlistenVal2D is 0, return TlistenVal without performing DP  21/08/05

    %#codegen
    coder.inline('never')
    
    t_axis = dt:dt:maxT;
    NSTEPS = length(t_axis); % Number of time-step grids
    SIGMA = sqrt(dt); % Standard deviation of the diffusion over one time step 

    NDV = (2*MAXDV/DVDELTA + 1); % Number of grids for the decision variable (DV)
    NBEL = NDV*NSTEPS;
    Nactions = 3;
    V_init = 0;
    rho = 0;
    
    policty_iter_flag = isempty(policy);
    
    uv = kappa * coherences;
    
    rewardSuccess = rewardPrm(1);
    rewardFailure = rewardPrm(2);
    rewardListen  = rewardPrm(3);
    rewardTimesUp = rewardPrm(4);
    
    % Mean inter-trial interval for valid trials (ITI_val) and canceled trials (ITI_inv)
    ITI_val = ITI(1);
    ITI_inv = ITI(2);
    
    %% preallocation
    Beliefs = zeros(NDV*NSTEPS,2);
    Beliefs2 = zeros(NDV*NSTEPS,2);
    
    qact = zeros(Nactions,1);
    V = zeros(NDV,NSTEPS);
    Vdelta = zeros(NDV,NSTEPS);
    Rewards = zeros(NBEL,3);
    updateOrder = zeros(NBEL,1);
    Q = zeros(NDV,NSTEPS,Nactions);
    bestQ = zeros(NDV,NSTEPS);
    
    % Assign some positive value to the state value at the origin
    % We will adjust the time cost (rho) until this state value becomes approximately zero.
    V((size(V,1)+1)/2,1) = 1;

    % Vector of the decison variable (dv)
    dv = -MAXDV:DVDELTA:MAXDV;

    cont = 0;
    for i = 1:NSTEPS
        for j = 1:NDV
            cont = cont+1;
            Beliefs(cont,1) = j; % dv
            Beliefs(cont,2) = i; % time
        end
    end
    cont = 0;
    for i = 1:NDV
        for j = 1:NSTEPS
            cont = cont+1;
            Beliefs2(cont,1) = i; % dv
            Beliefs2(cont,2) = j; % time
        end
    end
    
    % TlistenVal is a 3D transition matrix.
    % TlistenVal(i,j,k) is the transition probability of the decision variable (DV)
    % from DV(i) to DV(k) at time step t(j).
    % In order to pass it to a MEX function, it is saved as 2D matrix (TlistenVal2D).
    if isempty(TlistenVal2D)
        disp('creating TlistenVal');
        TlistenVal = makeTListenVal(NDV,NSTEPS,SIGMA,t_axis,dv,dt,uv); % same as the C code 07/17/14
    elseif numel(TlistenVal2D)==1 && ~all(TlistenVal2D(:))
        % if TlistenVal2D is 0, return TlistenVal without performing DP
        disp('creating TlistenVal and return');
        TlistenVal = makeTListenVal(NDV,NSTEPS,SIGMA,t_axis,dv,dt,uv); % same as the C code 07/17/14
        return
    else
        disp('using input TlistenVal');
        TlistenVal = reshape(TlistenVal2D,NDV,NSTEPS,[]);
    end

    rho1 = -1;
    rho2 = 1;

    while abs(V((size(V,1)+1)/2,1))>maxError

        rho = (rho1+rho2)/2;
        disp(rho);
        
        %% Initializations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for i=1:NBEL

            iDV = Beliefs(i,1);
            iT  = Beliefs(i,2); % time step

            % porbability that the given DV was generated by positive coherence
            p_plus  = calc_pplus(dv(iDV), iT, uv, SIGMA, t_axis);

            Rewards(i,1) = p_plus*rewardSuccess + (1-p_plus)*rewardFailure - ITI_val*rho; % reward for choosing plus (e.g. right target)
            Rewards(i,2) = (1-p_plus)*rewardSuccess + p_plus*rewardFailure  - ITI_val*rho; % reward for choosing minus (e.g., left target)
            Rewards(i,3) = rewardListen - rho*dt;

            if iT==NSTEPS
                Rewards(i,3) = rewardTimesUp - ITI_inv*rho;
            end
        end
        
        % Used as update order for the state values
        % Start updating for the last time step and then go backward in time (see below)
        cnt = 0;
        for j = NSTEPS:-1:1
            for i = 1:NBEL
                if Beliefs(i,2)==j
                    cnt = cnt+1;
                    updateOrder(cnt) = i;
                end
            end
        end
        
        %% main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        delta = 10;
        Vprev = zeros(NDV,NSTEPS);

        while delta > maxError
            
            % Because time flows in a single direction, the Bellman equations 
            % can be solved by backwards induction in a single pass.
            for iter = 1:NBEL
                iBelief = updateOrder(iter);
                iDV   = Beliefs(iBelief,1);
                iStep = Beliefs(iBelief,2);
                p_cancel = p_dl(iStep);
                % p_cancel become relevant and non-zero only for trials with provisional deadlines (Phase II).
                % p_cancel is the probability that the trail cancellation occurs at this time step, given it has not occured.
                % It is basically the harzard rate of trial cancellation.
                
                for iAction = 1:Nactions
                    qact(iAction) = Rewards(iBelief,iAction);
                    if ((iAction <= 2) || (iAction == 3 && iStep == NSTEPS)) % actions 0 to 3 are terminal
                        
                    elseif (iAction==3) % listen
                        temp = 0;
                        for j=1:NDV
                            p = TlistenVal(iDV,iStep,j);
                            temp = temp + p*(gamma*V(j,iStep+1)); % index for the next time step (istep+1)
                        end
                        qact(iAction) = qact(iAction) + (1-p_cancel)*temp + p_cancel*(-ITI_inv*rho);
                    end
                    Q(iDV,iStep,iAction) = qact(iAction); % Store action values (qact) for saving
                end

                [maxVal,maxInd] = max(qact);
                
                if policty_iter_flag
                    V(iDV,iStep) = maxVal;
                    bestQ(iDV,iStep) = maxInd;
                else
                    V(iDV,iStep) = qact(policy(iDV,iStep));
                    bestQ(iDV,iStep) = policy(iDV,iStep);
                end
            end

            Vdelta = abs(V - Vprev);
            Vprev  = V;  

            delta = max(max(Vdelta));

            V_init = V((size(V,1)+1)/2,1);
        end

        % Root finding, using the bisection method 
        % (Bertsekas, 1995; Drugowitsch et al., 2012).
        % We start by choosing two extreme values of ρ such that its true
        % value lies between them, and solve the Bellman equations using
        % value iteration, for the two extreme values and for their average
        % value. We calculate the difference between the reward rate we
        % assume and the one we get from deriving the optimal policy. We
        % then discard one of the two extreme values, the one for which the
        % error in the estimation of the reward rate has the same sign as
        % the central value. We repeat the process with the new interval
        % until the difference between the assumed reward rate and the one
        % obtained from solving the optimal policy is negligible.
        
        if V_init>0
            rho1 = rho;
        else
            rho2 = rho;
        end

    end
end

%%
function p_positive = calc_pplus(DV, timestep, uv, SIGMA, t_axis)
    NCOH = length(uv);
    p_coh = calc_pcoh(DV,timestep,uv, SIGMA, t_axis);

    % see positives
    p_positive = 0;
    for ci=1:NCOH
        if (uv(ci)>0)
            p_positive = p_positive + p_coh(ci);
        elseif (uv(ci)==0)
            p_positive = p_positive + p_coh(ci)/2;
        end
    end
end
%%
function p_coh = calc_pcoh(DV, timestep, uv, SIGMA, t_axis)
    log_p_coh = zeros(1,length(uv));
    ratio = zeros(1,length(uv));
    for ci=1:length(uv)
        dv_mean = uv(ci)*t_axis(timestep);
        dv_std =  sqrt(timestep)*SIGMA;
        log_p_coh(ci) = -1/2*((DV - dv_mean)/dv_std)^2;
    end
    m = ceil(length(uv)/2); % middle index
    for ci = 1:length(uv)
        ratio(ci) = exp(log_p_coh(ci) - log_p_coh(m)); % A/B
    end
    p_coh = ratio/sum(ratio);
end
%%
function TlistenVal = makeTListenVal(NDV,NSTEPS,SIGMA,t_axis,dv,dt,uv)
    NCOH = length(uv);
    TlistenVal = zeros(NDV,NSTEPS,NDV);
    for di=1:NDV
        for si=1:NSTEPS
            for k=1:NDV
                p_coh = calc_pcoh(dv(di), si, uv, SIGMA, t_axis);
                TlistenVal(di,si,k) = 0.0;
                for cohi=1:NCOH
                    TlistenVal(di,si,k) = TlistenVal(di,si,k) + p_coh(cohi)*normpdf(dv(k),dv(di)+uv(cohi)*dt,SIGMA);
                end
            end
            if sum(TlistenVal(di,si,:))>0
                TlistenVal(di,si,:) = TlistenVal(di,si,:)/sum(TlistenVal(di,si,:));
            else
                TlistenVal(di,si,:) = 1/NDV;
                disp([di,si]);
            end
        end
    end
end
