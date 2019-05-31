function [X,t]=SSA_Nk_channel_varying_v(t_fin, X0, Dt, Vpath)

%%% Code adapted from code for (Bruce, ABME 2009) available at
%%% http://www.ece.mcmaster.ca/~ibruce/ (authors website)

%%% INPUTS
% t_fin is final time of solution (ms)
% X0 is vector with NUMBER (so entries of X0 are integers) of channels in each state at beginning of
% the simulation (initial state of the channel), ordering of vector is as
% follows - X0=[n4 n3 n2 n1 n0]'
% Dt is the time step
% Vpath is the voltage path (so V varies over time) 

%%% OUTPUTS
% t is time
% X - X(:, i) is the NUMBER of channels in each state at time t(i) 
%(ordering of X(:, i) the same as X0).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise size of outputs
Nstep = ceil(t_fin/Dt); 
X=zeros(5, Nstep+1); %stores the voltage and proportion of open sodium channels at time intervals of length Dt.
t=zeros(Nstep+1,1);

% Vectors describing state transitions from current to next state
curr_state_k=[0 5 4 3 2 1 2 3 4];
next_state_k=[0 4 3 2 1 2 3 4 5];

% Set up initial conditions
X(:, 1)=X0;
y0=X0;%vector that tracks state of channel over time
t(1)=0;

% Main loop
for i=2:Nstep+1
    
    V=Vpath(i-1);% voltage constant over time step of length Dt
    
    % Transition rates, as function of the voltage
    alphan=(0.01*(10-V))/(exp((10-V)/10)-1);%original HH model
    betan=0.125*exp(-V/80);%original HH model
    
    zeta=[4*betan
        alphan+3*betan
        2*alphan+2*betan
        3*alphan+betan
        4*alphan];% Eq. (A3) of Mino et al. (ABME 2002)
    time=0;
    while time<Dt
              
        % Solve the system using SSA
        
        % Time till next reaction is an exponential random variable rate
        lambda=y0'*zeta;
        tau=exprnd((1/lambda));
        
        time=time+tau;
        
        if time>Dt
            break
        else
            % Deciding which reaction happens
            props_k=cumsum([0,...
                4*alphan*y0(5),...
                3*alphan*y0(4),...
                2*alphan*y0(3),...
                alphan*y0(2),...
                4*betan*y0(1),...
                3*betan*y0(2),...
                2*betan*y0(3),...
                betan*y0(4)]/lambda);
            
            % Determine which state transition has occurred
            ind_k = find(rand(1)<props_k,1);
            
            y0(curr_state_k(ind_k)) = y0(curr_state_k(ind_k)) - 1; % decrease number of channels in current state
            y0(next_state_k(ind_k)) = y0(next_state_k(ind_k)) + 1; % increase number of channels in next state
        end
    end
    
    % Update channel state (only do this at Dt time intervals)
    X(:, i)=y0;
    t(i)=t(i-1)+Dt;
    
end