function [X, t] = HH_RSDE_FE_sims (t_fin, X0, na0, Dt, I_amp, N, scaled)

%%% INPUTS
% t_fin is final time of solution (ms)
% X0 is the initial condition where X0(1) is initial voltage value and
% X0(2) is initial proportion of sodium channels in the open state
% na0 is vector with proportion of channels in each state at beginning of
% the simulation (initial state of the sodium channel), ordering of vector is as
% follows - X0=[m3h1 m2h1 m1h1 m0h1 m3h0 m2h0 m1h0 m0h0]'
% Dt is the time step
% I_amp is the amplitude of the stimulus current
% N is the number of channels
% Scaled indicates how number of channels is scaled (0 is constant membrane area, 1 is constant channel density)

%%% OUTPUTS
% t is time
% X - state of the system over time, X(1, :) is the voltage trace, X(2, :)
% is how proportion of sodium channels in open state varies over time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants for voltage equation taken from code available at http://www.ece.mcmaster.ca/~ibruce/ 

if scaled==0
    C=0.0714e-6;%membrane capacitance in uF [Cm = 0.0714nF] 
    R=1953.49e3;%membrane resistance in kOhms [Rm = 1.95349e3 MOhms]
else
    C = 0.0714e-6*N/1e3; % membrane capacitance in uF [Cm = 0.0714 pF]
    R = 1953.49e3/N*1e3; % membrane resistance in kOhms [Rm = 1.95349e3 MOhms]
end
g = 2.569e-8; %maximum sodium conductance in mS [gamma_Na = 25.69 pS/channel]
Ena=144;

% Initialise size of outputs
Nstep = ceil(t_fin/Dt);     %number of steps of size Dt=T/N
X=zeros(2, Nstep+1); 
na_channel=zeros(8, Nstep+1); %tracks proportion of sodium channels in each state
t=zeros(Nstep+1,1);

% Stimulus vector for monophasic pulse protocal
stim_duration = 0.1;   % millisecond (in membrane)
Istim = zeros(1,Nstep); % Stimulus vector
Istim(1:round(stim_duration/Dt-1)) = I_amp;

% Initial conditions
X(:, 1)=X0;
na_channel(:, 1)=na0;
t(1)=0;

% Brownian increments (scaled with the time step)
dW=sqrt(Dt)*randn(10, Nstep);

% Initialise matrices for the RSDE
E=[-1 0 0 0 0 0 0 0 0 -1
    1 -1 0 0 0 0 0 0 -1 0
    0 1 -1 0 0 0 0 -1 0 0
    0 0 1 0 0 0 -1 0 0 0
    0 0 0 -1 0 0 0 0 0 1
    0 0 0 1 -1 0 0 0 1 0
    0 0 0 0 1 -1 0 1 0 0];
E2=[E; 0 0 0 0 0 1 1 0 0 0];

% MAIN LOOP
for i=2:Nstep+1
    
    V=X(1, i-1); %voltage at current time step
    
    % Transition rates, as function of the voltage
    alpham=(1.872*(V-25.41))/(1-exp(-(V-25.410)/6.06));
    betam=(3.973*(21.001-V))/(1-exp((V-21.001)/9.41));
    alphah=(-0.549*(27.74+V))/(1-exp((V+27.74)/9.06));
    betah=22.57/(1+exp((56.0-V)/12.5));
    
    % Calculate drift term
    M=[-(betah+3*betam) alpham 0 0 alphah 0 0 0; 3*betam -(alpham+2*betam+betah) 2*alpham 0 0 alphah 0 0; 0 2*betam -(2*alpham+betam+betah) 3*alpham 0 0 alphah 0; 0 0 betam -(3*alpham+betah) 0 0 0 alphah];
    M=[M; betah 0 0 0 -(alphah+3*betam) alpham 0 0; 0 betah 0 0 3*betam -(alphah+alpham+2*betam) 2*alpham 0; 0 0 betah 0 0 2*betam -(2*alpham+betam+alphah) 3*alpham; 0 0 0 betah 0 0 betam -(alphah+3*alpham)];
    
    % Calculate part of the diffusion term
    S=sqrt([alpham*na_channel(2, i-1)+3*betam*na_channel(1, i-1) 2*alpham*na_channel(3, i-1)+2*betam*na_channel(2, i-1) 3*alpham*na_channel(4, i-1)+betam*na_channel(3, i-1) alpham*na_channel(6, i-1)+3*betam*na_channel(5, i-1) 2*alpham*na_channel(7, i-1)+2*betam*na_channel(6, i-1) 3*alpham*na_channel(8, i-1)+betam*na_channel(7, i-1) alphah*na_channel(4, i-1)+betah*na_channel(8, i-1) alphah*na_channel(3, i-1)+betah*na_channel(7, i-1) alphah*na_channel(2, i-1)+betah*na_channel(6, i-1) alphah*na_channel(1, i-1)+betah*na_channel(5, i-1)]);
    Noi=S'.*dW(:, i-1);
    
    % Calculate SDE without reflecting term, at next time step
    Xnext=na_channel(:, i-1)+M*na_channel(:, i-1)*Dt+(1/sqrt(N))*E2*Noi;
    
    % Check to see if point lies in correct domain
    if min(Xnext)>=0  && max(Xnext)<=1
        % Point lies in correct domain so set solution to be this value
        na_channel(:, i)=Xnext;
        
    else
        % Project X into the correct domain
        na_channel(:, i) = projsplx(Xnext);
        
    end
    
    % Calculate voltage at next time step using Euler method
    
    X(1, i)=X(1, i-1)+(1/C*((-X(1, i-1)/R)-g*(N*X(2, i-1))*(X(1, i-1)-Ena)+Istim(i-1)))*Dt;
    
    % Update
    X(2, i)=na_channel(1, i);
    t(i)=t(i-1)+Dt;
    
end