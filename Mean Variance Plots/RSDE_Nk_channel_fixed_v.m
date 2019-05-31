function [X, t] = RSDE_Nk_channel_fixed_v(t_fin, X0, Dt, v, N)

%%% INPUTS
% t_fin is final time of solution (ms)
% X0 is vector with proportion of channels in each state at beginning of
% the simulation (initial state of the channel), ordering of vector is as
% follows - X0=[n4 n3 n2 n1 n0]'
% Dt is the time step
% v is the voltage value, which remains fixed throughout the simulation
% N is the number of channels

%%% OUTPUTS
% t is time
% X - X(:, i) is the proportion of channels in each state at time t(i) 
%(ordering of X(:, i) the same as X0).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise size of outputs
Nstep = ceil(t_fin/Dt); % Number of steps of size Dt=T/N
X=zeros(5,Nstep+1); 
t=zeros(Nstep+1,1);

% Transition rates, as function of the voltage
alphan=(0.01*(10-v))/(exp((10-v)/10)-1);%original HH model
betan=0.125*exp(-v/80);%original HH model

% Initialise matrices for the RSDE
M=[-4*betan alphan 0 0 0
    4*betan -(alphan+3*betan) 2*alphan 0 0
    0 3*betan -(2*alphan+2*betan) 3*alphan 0
    0 0 2*betan -(3*alphan+betan) 4*alphan
    0 0 0 betan -4*alphan];
E=[1 0 0 0
    -1 1 0 0
    0 -1 1 0
    0 0 -1 1
    0 0 0 -1];

% Set up initial conditions
X(:, 1)=X0;
t(1)=0;

% Brownian increments (scaled with the time step)
dW = sqrt(Dt)*randn(4,Nstep);    

% Main loop
for i=2:Nstep+1
    
     % Calculate part of the diffusion term
     Noi=[sqrt((alphan*X(2, i-1)+4*betan*X(1, i-1)))*dW(1, i-1)
        sqrt((2*alphan*X(3, i-1)+3*betan*X(2, i-1)))*dW(2, i-1)
        sqrt((3*alphan*X(4, i-1)+2*betan*X(3, i-1)))*dW(3, i-1)
        sqrt((4*alphan*X(5, i-1)+betan*X(4, i-1)))*dW(4, i-1)];
    
    % Calculate SDE without reflecting term, at next time step
    Xnext=X(:, i-1)+(M*X(:, i-1))*Dt+(1/(sqrt(N)))*E*Noi;
    
    % Check to see if point lies in correct domain
    if min(Xnext)>=0  && max(Xnext)<=1
        
        % Point lies in correct domain so set solution to be this value
        X(:, i)=Xnext;
        
    else
        
        % Project X into the correct domain
        X(:, i) = projsplx(Xnext);
        
    end
    
   % Update the time vector       
   t(i)=t(i-1)+Dt;     
   
end