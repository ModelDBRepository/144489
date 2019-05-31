function [X, t] = HH_SDE_FE_sims (t_fin, X0, na0, Dt, I_amp, N, scaled)


%%% INPUTS
% t_fin is final time of solution (ms)
% X0 is the initial condition where X0(1) is initial voltage value and
% X0(2) is initial proportion of sodium channels in the open state
% na0 is vector with proportion of channels in each state at beginning of
% the simulation (initial state of the sodium channel), ordering of vector is as
% follows - X0=[m3h1 m2h1 m1h1 m0h1 m3h0 m2h0 m1h0]'
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
Nstep = ceil(t_fin/Dt);     % Number of steps of size Dt=T/N
X=zeros(2, Nstep+1); 
na_channel=zeros(7, Nstep+1); %tracks proportion of sodium channels in each state
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
dW=sqrt(Dt)*randn(7, Nstep);

% MAIN LOOP
for i=2:Nstep+1
    
    V=X(1, i-1); %voltage at current time step
    
    % Transition rates, as function of the voltage
    
    alpham=(1.872*(V-25.41))/(1-exp(-(V-25.410)/6.06));
    betam=(3.973*(21.001-V))/(1-exp((V-21.001)/9.41));
    alphah=(-0.549*(27.74+V))/(1-exp((V+27.74)/9.06));
    betah=22.57/(1+exp((56.0-V)/12.5));
    
    % Calculate drift term
    Mna=[-(betah+3*betam) alpham 0 0 alphah 0 0; 3*betam -(alpham+2*betam+betah) 2*alpham 0 0 alphah 0; 0 2*betam -(2*alpham+betam+betah) 3*alpham 0 0 alphah; -alphah -alphah betam-alphah -(3*alpham+betah+alphah) -alphah -alphah -alphah];
    Mna=[Mna; betah 0 0 0 -(alphah+3*betam) alpham 0; 0 betah 0 0 3*betam -(alphah+alpham+2*betam) 2*alpham; -3*alpham -3*alpham betah-3*alpham -3*alpham -3*alpham 2*betam-3*alpham -(2*alpham+betam+alphah+3*alpham)];
    con_na=[ 0;  0;  0;  alphah; 0;  0; 3*alpham];
    
    % Calculate steady state of the channel, used to calulate diffusion term
    na_channel_st=1/((alpham+betam)^3*(alphah+betah))*[alpham^3*alphah*nchoosek(3,3)
        alpham^2*betam*alphah*nchoosek(3,2)
        alpham*betam^2*alphah*nchoosek(3,1)
        betam^3*alphah*nchoosek(3,0)
        alpham^3*betah*nchoosek(3,3)
        alpham^2*betam*betah*nchoosek(3,2)
        alpham*betam^2*betah*nchoosek(3,1)
        betam^3*betah*nchoosek(3,0)];
    
    % Calculate part of the diffusion term
    Dna=(1/N)*[alphah*na_channel_st(5)+alpham*na_channel_st(2)+(3*betam+betah)*na_channel_st(1) -(alpham*na_channel_st(2)+3*betam*na_channel_st(1)) 0 0 -(alphah*na_channel_st(5)+betah*na_channel_st(1)) 0 0
        -(alpham*na_channel_st(2)+3*betam*na_channel_st(1)) alphah*na_channel_st(6)+2*alpham*na_channel_st(3)+(alpham+2*betam+betah)*na_channel_st(2)+3*betam*na_channel_st(1) -2*(alpham*na_channel_st(3)+betam*na_channel_st(2)) 0 0 -(alphah*na_channel_st(6)+betah*na_channel_st(2)) 0
        0 -2*(alpham*na_channel_st(3)+betam*na_channel_st(2)) alphah*na_channel_st(7)+3*alpham*na_channel_st(4)+(2*alpham+betam+betah)*na_channel_st(3)+2*betam*na_channel_st(2) -(3*alpham*na_channel_st(4)+betam*na_channel_st(3)) 0 0 -(alphah*na_channel_st(7)+betah*na_channel_st(3))
        0 0 -(3*alpham*na_channel_st(4)+betam*na_channel_st(3)) alphah*na_channel_st(8)+(3*alpham+betah)*na_channel_st(4)+betam*na_channel_st(3) 0 0 0
        -(alphah*na_channel_st(5)+betah*na_channel_st(1)) 0 0 0 alpham*na_channel_st(6)+(3*betam+alphah)*na_channel_st(5)+betah*na_channel_st(1) -(alpham*na_channel_st(6)+3*betam*na_channel_st(5)) 0
        0 -(alphah*na_channel_st(6)+betah*na_channel_st(2)) 0 0 -(alpham*na_channel_st(6)+3*betam*na_channel_st(5)) 2*alpham*na_channel_st(7)+(alpham+2*betam+alphah)*na_channel_st(6)+3*betam*na_channel_st(5)+betah*na_channel_st(2) -2*(alpham*na_channel_st(7)+betam*na_channel_st(6))
        0 0 -(alphah*na_channel_st(7)+betah*na_channel_st(3)) 0 0 -2*(alpham*na_channel_st(7)+betam*na_channel_st(6)) 3*alpham*na_channel_st(8)+(2*alpham+betam+alphah)*na_channel_st(7)+2*betam*na_channel_st(6)+betah*na_channel_st(3)];
    
    % Calculate SDE at next time step using Euler Maruyama method
    na_channel(:, i)=na_channel(:, i-1)+(Mna*na_channel(:, i-1)+con_na)*Dt+sqrtm(Dna)*dW(:, i-1);
    
    % Calculate voltage at next time step using Euler method
    X(1, i)=X(1, i-1)+(1/C*((-X(1, i-1)/R)-g*(N*X(2, i-1))*(X(1, i-1)-Ena)+Istim(i-1)))*Dt;
    
    % Update
    X(2, i)=na_channel(1, i);
    t(i)=t(i-1)+Dt;
    
end