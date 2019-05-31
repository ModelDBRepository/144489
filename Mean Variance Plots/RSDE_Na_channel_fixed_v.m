function [X, t] = RSDE_Na_channel_fixed_v(t_fin, X0, Dt, v, N)

%%% INPUTS
% t_fin is final time of solution (ms)
% X0 is vector with proportion of channels in each state at beginning of
% the simulation (initial state of the channel), ordering of vector is as
% follows - X0=[m3h1 m2h1 m1h1 m0h1 m3h0 m2h0 m1h0 m0h0]'
% Dt is the time step
% v is the voltage value, which remains fixed throughout the simulation
% N is the number of channels

%%% OUTPUTS
% t is time
% X - X(:, i) is the proportion of channels in each state at time t(i) 
%(ordering of X(:, i) the same as X0).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise size of outputs
Nstep = ceil(t_fin/Dt);     % Number of steps of size Dt=T/N
X=zeros(8,Nstep+1);
t=zeros(Nstep+1,1);

% Transition rates, as function of the voltage
alpham=0.1*((25-v)/(exp((25-v)/10)-1));
betam=4*exp(-v/18);
alphah=0.07*exp(-v/20);
betah=1/(exp((30-v)/10)+1);

% Initialise matrices for the RSDE
M=[-(betah+3*betam) alpham 0 0 alphah 0 0 0; 3*betam -(alpham+2*betam+betah) 2*alpham 0 0 alphah 0 0; 0 2*betam -(2*alpham+betam+betah) 3*alpham 0 0 alphah 0; 0 0 betam -(3*alpham+betah) 0 0 0 alphah];
M=[M; betah 0 0 0 -(alphah+3*betam) alpham 0 0; 0 betah 0 0 3*betam -(alphah+alpham+2*betam) 2*alpham 0; 0 0 betah 0 0 2*betam -(2*alpham+betam+alphah) 3*alpham; 0 0 0 betah 0 0 betam -(alphah+3*alpham)];
E=[-1 0 0 0 0 0 0 0 0 -1
    1 -1 0 0 0 0 0 0 -1 0
    0 1 -1 0 0 0 0 -1 0 0
    0 0 1 0 0 0 -1 0 0 0
    0 0 0 -1 0 0 0 0 0 1
    0 0 0 1 -1 0 0 0 1 0
    0 0 0 0 1 -1 0 1 0 0
    0 0 0 0 0 1 1 0 0 0];

% Set up initial conditions
X(:, 1)=X0;
t(1)=0;

% Brownian increments (scaled with the time step)
dW = sqrt(Dt)*randn(10,Nstep);    

% Main loop
for i=2:Nstep+1
    
    % Calculate part of the diffusion term
    S=sqrt([alpham*X(2, i-1)+3*betam*X(1, i-1) 2*alpham*X(3, i-1)+2*betam*X(2, i-1) 3*alpham*X(4, i-1)+betam*X(3, i-1) alpham*X(6, i-1)+3*betam*X(5, i-1) 2*alpham*X(7, i-1)+2*betam*X(6, i-1) 3*alpham*X(8, i-1)+betam*X(7, i-1) betah*X(4, i-1)+alphah*X(8, i-1) betah*X(3, i-1)+alphah*X(7, i-1) betah*X(2, i-1)+alphah*X(6, i-1) betah*X(1, i-1)+alphah*X(5, i-1)]);
    Noi=S'.*dW(:, i-1);
    
    % Calculate SDE without reflecting term, at next time step
    Xnext=X(:, i-1)+M*X(:, i-1)*Dt+(1/sqrt(N))*E*Noi;
    
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
