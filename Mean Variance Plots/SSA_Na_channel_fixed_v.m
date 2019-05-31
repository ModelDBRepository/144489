function [X,t]=SSA_Na_channel_fixed_v(t_fin, X0, Dt, v)

%%% Code adapted from code for (Bruce, ABME 2009) available at
%%% http://www.ece.mcmaster.ca/~ibruce/ (authors website)

%%% INPUTS
% t_fin is final time of solution (ms)
% X0 is vector with NUMBER (so entries of X0 are integers) of channels in each state at beginning of
% the simulation (initial state of the channel), ordering of vector is as
% follows - X0=[m3h1 m2h1 m1h1 m0h1 m3h0 m2h0 m1h0 m0h0]'
% Dt is the time step
% v is the voltage value, which remains fixed throughout the simulation

%%% OUTPUTS
% t is time
% X - X(:, i) is the NUMBER of channels in each state at time t(i) 
%(ordering of X(:, i) the same as X0).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise size of outputs
Nstep = ceil(t_fin/Dt); 
X=zeros(8, Nstep+1); %stores the voltage and proportion of open sodium channels at time intervals of length Dt.
t=zeros(Nstep+1,1);

% Transition rates, as function of the voltage
alpham=0.1*((25-v)/(exp((25-v)/10)-1));%original HH model
betam=4*exp(-v/18);%original HH model
alphah=0.07*exp(-v/20);%original HH model
betah=1/(exp((30-v)/10)+1);%original HH model


% Set up initial conditions
X(:, 1)=X0;
mh=X(:, 1); %vector that tracks state of channel over time
t(1)=0;

zeta=[ 3*betam + betah
       2*betam + betah + alpham
       betam + betah + 2*alpham
       3*alpham + betah
       3*betam + alphah
       2*betam + alphah + alpham
       betam + alphah + 2*alpham
       3*alpham + alphah]; % Eq. (A3) of Mino et al. (ABME 2002)
   
% Vectors describing state transitions from current to next state
curr_state=[0 4 3 2 8 7 6 1 2 3 5 6 7 8 7 6 5 4 3 2 1];
next_state=[0 3 2 1 7 6 5 2 3 4 6 7 8 4 3 2 1 8 7 6 5];    

% Main loop    
for i=2:Nstep+1
    
    time=0;
    
    while time<Dt
               
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
    
    % Update channel state (only do this at Dt time intervals)
    X(:, i)=mh;
    t(i)=t(i-1)+Dt;
    
end