function [X, t] = HH_SSA_FE_sims (t_fin, X0, X0_na, Dt, I_amp, scaled, Nna)

%%% Code adapted from code for (Bruce, ABME 2009) available at
%%% http://www.ece.mcmaster.ca/~ibruce/ (authors website)

%%% INPUTS
% t_fin is final time of solution (ms)
% X0 is the initial condition where X0(1) is initial voltage value and
% X0(2) is initial number of sodium channels in the open state
% X0na is vector with NUMBER (so entries of X0 are integers) of channels in each state at beginning of
% the simulation (initial state of the sodium channel), ordering of vector is as
% follows - X0=[m3h1 m2h1 m1h1 m0h1 m3h0 m2h0 m1h0 m0h0]'
% Dt is the time step
% I_amp is the amplitude of the stimulus current
% Scaled indicates how number of channels is scaled (0 is constant membrane area, 1 is constant channel density)
% Nna is the number of channels


%%% OUTPUTS
% t is time
% X - state of the system over time, X(1, :) is the voltage trace, X(2, :)
% is how number of sodium channels in open state varies over time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants for voltage equation taken from code available at http://www.ece.mcmaster.ca/~ibruce/ 

if scaled==0
    C=0.0714e-6;%membrane capacitance in uF [Cm = 0.0714nF] 
    R=1953.49e3;%membrane resistance in kOhms [Rm = 1.95349e3 MOhms] 
else
    C = 0.0714e-6*Nna/1e3; % membrane capacitance in uF [Cm = 0.0714 pF]
    R = 1953.49e3/Nna*1e3; % membrane resistance in kOhms [Rm = 1.95349e3 MOhms]
end
g = 2.569e-8; %maximum sodium conductance in mS [gamma_Na = 25.69 pS/channel]
Ena=144;

% Initialise size of outputs
Nstep = ceil(t_fin/Dt); 
X=zeros(2, Nstep+1); % Stores the voltage and proportion of open sodium channels at time intervals of length Dt.
t=zeros(Nstep+1,1);

% Stimulus vector for monophasic pulse protocal
stim_duration = 0.1;   % millisecond (in membrane)
Istim = zeros(1,Nstep); % Stimulus vector
Istim(1:round(stim_duration/Dt-1)) = I_amp;

% Set up initial conditions
X(:, 1)=X0;
mh=X0_na; %vector that tracks state of channel over time
t(1)=0;

% MAIN LOOP
for i=2:Nstep+1
    
    V=X(1, i-1); %voltage at current time step
    
    % Transition rates, as function of the voltage
    alpham=(1.872*(V-25.41))/(1-exp(-(V-25.410)/6.06));
    betam=(3.973*(21.001-V))/(1-exp((V-21.001)/9.41));
    alphah=(-0.549*(27.74+V))/(1-exp((V+27.74)/9.06));
    betah=22.57/(1+exp((56.0-V)/12.5));
    
    % Calculate transition rates of "escapes" from state m_ih_j (i=0.1.2.3; j=0,1)
    zeta=[ 3*betam + betah
        2*betam + betah + alpham
        betam + betah + 2*alpham
        3*alpham + betah
        3*betam + alphah
        2*betam + alphah + alpham
        betam + alphah + 2*alpham
        3*alpham + alphah]; % Eq. (A3) of Mino et al. (ABME 2002)
    
    time=0;
    while time<Dt
        % Vectors describing state transitions from current to next state
        curr_state=[0 4 3 2 8 7 6 1 2 3 5 6 7 8 7 6 5 4 3 2 1];
        next_state=[0 3 2 1 7 6 5 2 3 4 6 7 8 4 3 2 1 8 7 6 5];
        
        % Solve the system using SSA
        
        % Time till next reaction is an exponential random variable rate
        lambda=mh'*zeta;
        tau=exprnd((1/lambda));
        time=time+tau;
        
        if time>Dt
            break
        else
            % Deciding which reaction happens            
            props=cumsum([0,...
                3*alpham*mh(4),...
                2*alpham*mh(3),...
                alpham*mh(2),...
                3*alpham*mh(8),...
                2*alpham*mh(7),...
                alpham*mh(6),...
                3*betam*mh(1),...
                2*betam*mh(2),...
                betam*mh(3),...
                3*betam*mh(5),...
                2*betam*mh(6),...
                betam*mh(7),...
                alphah*mh(8),...
                alphah*mh(7),...
                alphah*mh(6),...
                alphah*mh(5),...
                betah*mh(4),...
                betah*mh(3),...
                betah*mh(2),...
                betah*mh(1)]/lambda);
            
            % Determine which state transition has occurred
            ind = find(rand(1)<props,1);
            
            mh(curr_state(ind)) = mh(curr_state(ind)) - 1; % decrease number of channels in current state
            mh(next_state(ind)) = mh(next_state(ind)) + 1; % increase number of channels in next state
            
        end
        
    end
    
    % Calculate voltage at next time step using Euler method
    
    X(1, i)=V+(1/C*((-V/R)-g*X(2, i-1)*(V-Ena)+Istim(i-1)))*Dt;
    
    % Update values in X and t
    
    X(2, i)=mh(1);
    t(i)=t(i-1)+Dt;
    
end