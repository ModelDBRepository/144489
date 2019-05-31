function [X, t] = RSDE_Nk_channel_varying_v(t_fin, X0, Dt, N, Vpath)

%%% INPUTS
% t_fin is final time of solution (ms)
% X0 is vector with proportion of channels in each state at beginning of
% the simulation (initial state of the channel), ordering of vector is as
% follows - X0=[n4 n3 n2 n1 n0]'
% Dt is the time step
% N is the number of channels
% Vpath is the voltage path (so V varies over time)

%%% OUTPUTS
% t is time
% X - X(:, i) is the proportion of channels in each state at time t(i) 
%(ordering of X(:, i) the same as X0).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise size of outputs
Nstep = ceil(t_fin/Dt);     % Number of steps of size Dt=T/N
X=zeros(5,Nstep+1); %solution for orignal langevin
t=zeros(Nstep+1,1);

% Initialise matrices for the RSDE
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
    
    % Set voltage to value at current time step
    V=Vpath(i);
    
    % Calculate transition rates at current timestep, as function of the voltag
    an=(0.01*(10-V))/(exp((10-V)/10)-1);%original HH model
    bn=0.125*exp(-V/80);%original HH model
   
    % Calculate drift term
    M=[-4*bn an 0 0 0
        4*bn -(an+3*bn) 2*an 0 0
        0 3*bn -(2*an+2*bn) 3*an 0
        0 0 2*bn -(3*an+bn) 4*an
        0 0 0 bn -4*an];
     
    % Calculate diffusion term
    Noi=[sqrt((an*X(2, i-1)+4*bn*X(1, i-1)))*dW(1, i-1)
        sqrt((2*an*X(3, i-1)+3*bn*X(2, i-1)))*dW(2, i-1)
        sqrt((3*an*X(4, i-1)+2*bn*X(3, i-1)))*dW(3, i-1)
        sqrt((4*an*X(5, i-1)+bn*X(4, i-1)))*dW(4, i-1)];
    
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