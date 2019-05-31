function [jitter_sde_d, meanlt_sde_d, jitter_ssa_d, meanlt_ssa_d, jitter_rsde_d, meanlt_rsde_d, jitter_sde_a, meanlt_sde_a, jitter_ssa_a, meanlt_ssa_a, jitter_rsde_a, meanlt_rsde_a] = AP_statistics_sims ();

%%% Code adapted from code for (Bruce, ABME 2009) available at
%%% http://www.ece.mcmaster.ca/~ibruce/ (authors website)

%%% OUTPUTS
% Figures 9 and 10 in Dangerfield et al. Phys Rev E, Vol. 85 Issue 5 2012

clc
clear all
close all

Nna=[100 500 1000 5000 10000]'; %number of sodium channels
Dt=0.001; % time step used in simulation
t_fin=2; % length of simulation
%N_sim=1000; % number of simulations
N_sim=100;

% Run ODE describing Markov state of channel to steady state for V=0 to
% obtain initial condition for Na channel in RSDE and SDE sims
v=0;
x0=[1/8; 1/8; 1/8; 1/8; 1/8; 1/8; 1/8; 1/8];
tspan=[0, 100];
[to_Na, xo_Na]=ode45(@(t, x) mh_gate_HH_ODE(t, x, v),tspan,x0);
ic=xo_Na(end, :);

% Initial condition for RSDE and SDE sims
X0_sde=[0; ic(1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                   Constant Membrane Area                        %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scaled=0;

% tSimulus amplitude - values taken from findIthRS.m available at http://www.ece.mcmaster.ca/~ibruce/
I_amp=zeros(length(Nna), 17);
I_amp(1, :)=linspace(28,40,17); % For Na_max = 100, R & C not scaled
I_amp(2, :)=linspace(20,28,17); % For Na_max = 500, R & C not scaled
I_amp(3, :)=linspace(16,26,17); % For Na_max = 1,000
I_amp(4, :)=linspace(14,19,17); % For Na_max = 5,000, R & C not scaled
I_amp(5, :)=linspace(11,17,17); % For Na_max = 10,000, R & C not scaled
I_amp=I_amp*1e-6;

%%% Rows correspond to numbers of channels, columns to different stim amp

% Firing efficiencies
FE_ssa_a=zeros(length(Nna), length(I_amp(1, :))); % calculated using SSA
FE_rsde_a=zeros(length(Nna), length(I_amp(1, :))); % calculated using RSDE
FE_sde_a=zeros(length(Nna), length(I_amp(1, :))); % calculated using SDE with equilibrium noise term
% Mean latency
meanlt_ssa_a=zeros(length(Nna), length(I_amp)); % calculated using SSA
meanlt_rsde_a=zeros(length(Nna), length(I_amp)); % calculated using RSDE
meanlt_sde_a=zeros(length(Nna), length(I_amp)); % calculated using SDE with equilibrium noise term
% Jitter
jitter_ssa_a=zeros(length(Nna), length(I_amp)); % calculated using SSA
jitter_rsde_a=zeros(length(Nna), length(I_amp)); % calculated using RSDE
jitter_sde_a=zeros(length(Nna), length(I_amp)); % calculated using SDE with equilibrium noise term

for k=1:length(Nna)
    % Row is simulation current, column is simulation number
    % Store firing occurances
    fired_ssa_a=zeros(length(I_amp), N_sim);
    fired_rsde_a=zeros(length(I_amp), N_sim);
    fired_sde_a=zeros(length(I_amp), N_sim);
    % Latency
    late_ssa_a=zeros(length(I_amp), N_sim);
    late_rsde_a=zeros(length(I_amp), N_sim);
    late_sde_a=zeros(length(I_amp), N_sim);
    
    % Initial condition for SSA simulations
    X0_na_in=round(Nna(k)*ic(1:7))';
    X0_na=[X0_na_in; Nna(k)-sum(X0_na_in)];
    X0=[0; X0_na(1)];
    
    for i=1:length(I_amp)
        i
        k
        for j=1:N_sim
            [X_ssa, t_ssa] = HH_SSA_FE_sims (t_fin, X0, X0_na, Dt, I_amp(k, i), scaled, Nna(k));
            [X_rsde, t_rsde] = HH_RSDE_FE_sims (t_fin, X0_sde, ic', Dt, I_amp(k, i), Nna(k), scaled);
            [X_sde, t_sde] = HH_SDE_FE_sims (t_fin, X0_sde, (ic(1:7))', Dt, I_amp(k,i), Nna(k), scaled);
            
            % Determine in each case whether or not spike has occured
            
            [diff_volt_1, pos_1]=max(X_ssa(1,:));
            [diff_volt_2, pos_2]=max(X_rsde(1,:));
            [diff_volt_3, pos_3]=max(X_sde(1,:));
            
            if diff_volt_1>=45
                fired_ssa_a(i, j)=1;
                late_ssa_a(i,j)=t_ssa(pos_1);
            else
                late_ssa_a(i,j)=NaN;
            end
            
            if diff_volt_2>=45
                fired_rsde_a(i, j)=1;
                late_rsde_a(i,j)=t_rsde(pos_2);
            else
                late_rsde_a(i,j)=NaN;
            end
            
            if diff_volt_3>=45
                fired_sde_a(i, j)=1;
                late_sde_a(i,j)=t_sde(pos_3);
            else
                late_sde_a(i,j)=NaN;
            end
        end
        
        FE_ssa_a(k,i)=sum(fired_ssa_a(i, :))/N_sim;
        FE_rsde_a(k,i)=sum(fired_rsde_a(i, :))/N_sim;
        FE_sde_a(k,i)=sum(fired_sde_a(i, :))/N_sim;
        
        if FE_ssa_a(k,i)>0
            meanlt_ssa_a(k,i) = mean(late_ssa_a(i, find(~isnan(late_ssa_a(i, :)))));
            jitter_ssa_a(k,i) = std(late_ssa_a(i, find(~isnan(late_ssa_a(i, :)))));
        else
            meanlt_ssa_a(k,i) = NaN;
            jitter_ssa_a(k,i) = NaN;
        end
        
        if FE_rsde_a(k,i)>0
            meanlt_rsde_a(k,i) = mean(late_rsde_a(i, find(~isnan(late_rsde_a(i, :)))));
            jitter_rsde_a(k,i) = std(late_rsde_a(i, find(~isnan(late_rsde_a(i, :)))));
        else
            meanlt_rsde_a(k,i) = NaN;
            jitter_rsde_a(k,i) = NaN;
        end
        
        if FE_sde_a(k,i)>0
            meanlt_sde_a(k,i) = mean(late_sde_a(i, find(~isnan(late_sde_a(i, :)))));
            jitter_sde_a(k,i) = std(late_sde_a(i, find(~isnan(late_sde_a(i, :)))));
        else
            meanlt_sde_a(k,i) = NaN;
            jitter_sde_a(k,i) = NaN;
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                   Calculating Ith and RS                        %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ith_ssa_a=zeros(length(Nna), 1);
Ith_rsde_a=zeros(length(Nna), 1);
Ith_sde_a=zeros(length(Nna), 1);
RS_ssa_a=zeros(length(Nna), 1);
RS_rsde_a=zeros(length(Nna), 1);
RS_sde_a=zeros(length(Nna), 1);
I_amp=I_amp*1e6;
for k=1:length(Nna)
    
    Yssa = fminsearch(@(x) intgausserr(x,I_amp(k, :),FE_ssa_a(k, :)),[mean(I_amp(k, :)) 0.03*mean(I_amp(k, :))]);
    Ith_ssa_a(k) = Yssa(1);
    RS_ssa_a(k) = Yssa(2)/Yssa(1); % Eq. (17) from Bruce 2009
    
    Yrsde = fminsearch(@(x) intgausserr(x,I_amp(k, :),FE_rsde_a(k, :)),[mean(I_amp(k, :)) 0.03*mean(I_amp(k, :))]);
    Ith_rsde_a(k) = Yrsde(1);
    RS_rsde_a(k) = Yrsde(2)/Yrsde(1); % Eq. (17) from Bruce 2009
    
    Ysde = fminsearch(@(x) intgausserr(x,I_amp(k, :),FE_sde_a(k, :)),[mean(I_amp(k, :)) 0.03*mean(I_amp(k, :))]);
    Ith_sde_a(k) = Ysde(1);
    RS_sde_a(k) = Ysde(2)/Ysde(1); % Eq. (17) from Bruce 2009
    
end

clear 'Iamp'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                   Constant Channel Density                      %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scaled=1;

% Stimulus amplitude - values taken from findIthRS.m available at http://www.ece.mcmaster.ca/~ibruce/
I_amp=zeros(length(Nna), 17);
I_amp(1, :)=linspace(12,30,17)*(100)/1e3; % For Na_max = 100, scaled R & C
I_amp(2, :)=linspace(18,26,17)*(500)/1e3; % For Na_max = 500, scaled R & C
I_amp(3, :)=linspace(16,26,17)*(1000)/1e3; % For Na_max = 1,000
I_amp(4, :)=linspace(20,23,17)*(5000)/1e3; % For Na_max = 5,000, scaled R & C
I_amp(5, :)=linspace(20.5,23,17)*(10000)/1e3; % For Na_max = 10,000, scaled R & C
I_amp=I_amp*1e-6;

%%% Rows correspond to numbers of channels, columns to different stim amp

% Firing efficiencies
FE_ssa_d=zeros(length(Nna), length(I_amp(1, :)));
FE_rsde_d=zeros(length(Nna), length(I_amp(1, :)));
FE_sde_d=zeros(length(Nna), length(I_amp(1, :)));
% Mean Latency
meanlt_ssa_d=zeros(length(Nna), length(I_amp(1, :)));
meanlt_rsde_d=zeros(length(Nna), length(I_amp(1, :)));
meanlt_sde_d=zeros(length(Nna), length(I_amp(1, :)));
% Jitter
jitter_ssa_d=zeros(length(Nna), length(I_amp(1, :)));
jitter_rsde_d=zeros(length(Nna), length(I_amp(1, :)));
jitter_sde_d=zeros(length(Nna), length(I_amp(1, :)));

for k=1:length(Nna)
    % Row is simulation current, column is simulation number
    % Store firing occurances
    fired_ssa_d=zeros(length(I_amp), N_sim);
    fired_rsde_d=zeros(length(I_amp), N_sim);
    fired_sde_d=zeros(length(I_amp), N_sim);
    % Latency
    late_ssa_d=zeros(length(I_amp), N_sim);
    late_rsde_d=zeros(length(I_amp), N_sim);
    late_sde_d=zeros(length(I_amp), N_sim);
    
    % Initial condition for SSA simulations
    X0_na_in=round(Nna(k)*ic(1:7))';
    X0_na=[X0_na_in; Nna(k)-sum(X0_na_in)];
    X0=[0; X0_na(1)];
    for i=1:length(I_amp)
        i
        k
        for j=1:N_sim
           
            [X_ssa, t_ssa] = HH_SSA_FE_sims (t_fin, X0, X0_na, Dt, I_amp(k, i), scaled, Nna(k));
            [X_rsde, t_rsde] = HH_RSDE_FE_sims (t_fin, X0_sde, ic', Dt, I_amp(k, i), Nna(k), scaled);
            [X_sde, t_sde] = HH_SDE_FE_sims (t_fin, X0_sde, (ic(1:7))', Dt, I_amp(k,i), Nna(k), scaled);
            
            % Determine in each case whether or not spike has occured
            
            [diff_volt_1, pos_1]=max(X_ssa(1,:));
            [diff_volt_2, pos_2]=max(X_rsde(1,:));
            [diff_volt_3, pos_3]=max(X_sde(1,:));
            
            if diff_volt_1>=45
                fired_ssa_d(i, j)=1;
                late_ssa_d(i,j)=t_ssa(pos_1);
            else
                late_ssa_d(i,j)=NaN;
            end
            
            if diff_volt_2>=45
                fired_rsde_d(i, j)=1;
                late_rsde_d(i,j)=t_rsde(pos_2);
            else
                late_rsde_d(i,j)=NaN;
            end
            
            if diff_volt_3>=45
                fired_sde_d(i, j)=1;
                late_sde_d(i,j)=t_sde(pos_3);
            else
                late_sde_d(i,j)=NaN;
            end
            
            FE_ssa_d(k,i)=sum(fired_ssa_d(i, :))/N_sim;
            FE_rsde_d(k,i)=sum(fired_rsde_d(i, :))/N_sim;
            FE_sde_d(k,i)=sum(fired_sde_d(i, :))/N_sim;
            
            if FE_ssa_d(k,i)>0
                meanlt_ssa_d(k,i) = mean(late_ssa_d(i, find(~isnan(late_ssa_d(i, :)))));
                jitter_ssa_d(k,i) = std(late_ssa_d(i, find(~isnan(late_ssa_d(i, :)))));
            else
                meanlt_ssa_d(k,i) = NaN;
                jitter_ssa_d(k,i) = NaN;
            end
            
            if FE_rsde_d(k,i)>0
                meanlt_rsde_d(k,i) = mean(late_rsde_d(i, find(~isnan(late_rsde_d(i, :)))));
                jitter_rsde_d(k,i) = std(late_rsde_d(i, find(~isnan(late_rsde_d(i, :)))));
            else
                meanlt_rsde_d(k,i) = NaN;
                jitter_rsde_d(k,i) = NaN;
            end
            
            if FE_sde_d(k,i)>0
                meanlt_sde_d(k,i) = mean(late_sde_d(i, find(~isnan(late_sde_d(i, :)))));
                jitter_sde_d(k,i) = std(late_sde_d(i, find(~isnan(late_sde_d(i, :)))));
            else
                meanlt_sde_d(k,i) = NaN;
                jitter_sde_d(k,i) = NaN;
            end
            
            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                   Calculating Ith and RS                        %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ith_ssa_d=zeros(length(Nna), 1);
Ith_rsde_d=zeros(length(Nna), 1);
Ith_sde_d=zeros(length(Nna), 1);
RS_ssa_d=zeros(length(Nna), 1);
RS_rsde_d=zeros(length(Nna), 1);
RS_sde_d=zeros(length(Nna), 1);
I_amp=I_amp*1e6;
for k=1:length(Nna)
    
    Yssa = fminsearch(@(x) intgausserr(x,I_amp(k, :),FE_ssa_d(k, :)),[mean(I_amp(k, :)) 0.03*mean(I_amp(k, :))]);
    Ith_ssa_d(k) = Yssa(1);
    RS_ssa_d(k) = Yssa(2)/Yssa(1); % Eq. (17) from Bruce 2009
    
    Yrsde = fminsearch(@(x) intgausserr(x,I_amp(k, :),FE_rsde_d(k, :)),[mean(I_amp(k, :)) 0.03*mean(I_amp(k, :))]);
    Ith_rsde_d(k) = Yrsde(1);
    RS_rsde_d(k) = Yrsde(2)/Yrsde(1); % Eq. (17) from Bruce 2009
    
    Ysde= fminsearch(@(x) intgausserr(x,I_amp(k, :),FE_sde_d(k, :)),[mean(I_amp(k, :)) 0.03*mean(I_amp(k, :))]);
    Ith_sde_d(k) = Ysde(1);
    RS_sde_d(k) = Ysde(2)/Ysde(1); % Eq. (17) from Bruce 2009
    
end

clear 'Iamp'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                        FIGURES                                  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fig. 9 in Dangerfield et al. Phys Rev E 2012

figure
subplot(2, 2, 1)
semilogx(Nna, Ith_ssa_d, 'r', 'LineWidth', 3,'MarkerSize', 5)
hold on
semilogx(Nna, Ith_rsde_d, 'o--b', 'LineWidth', 3,'MarkerSize', 15)
semilogx(Nna, Ith_sde_d, '^:k', 'LineWidth', 3,'MarkerSize', 5)
hold off
xlabel 'N_{Na}'
ylabel 'Ith'
title ('constant density')
subplot(2, 2, 3)
semilogx(Nna, RS_ssa_d, 'r', 'LineWidth', 3,'MarkerSize', 5)
hold on
semilogx(Nna, RS_rsde_d, 'o--b', 'LineWidth', 3,'MarkerSize', 15)
semilogx(Nna, RS_sde_d, '^:k', 'LineWidth', 3,'MarkerSize', 5)
hold off
xlabel 'N_{Na}'
ylabel 'RS'
subplot(2, 2, 2)
semilogx(Nna, Ith_ssa_a, 'r','LineWidth', 3,'MarkerSize', 5)
hold on
semilogx(Nna, Ith_rsde_a, 'o--b','LineWidth', 3,'MarkerSize', 15)
semilogx(Nna, Ith_sde_a, '^:k','LineWidth', 3,'MarkerSize', 5)
hold off
xlabel 'N_{Na}'
ylabel 'Ith'
title ('constant area')
subplot(2, 2, 4)
semilogx(Nna, RS_ssa_a, 'r','LineWidth', 3,'MarkerSize', 5)
hold on
semilogx(Nna, RS_rsde_a, 'o--b','LineWidth', 3,'MarkerSize', 15)
semilogx(Nna, RS_sde_a, '^:k','LineWidth', 3,'MarkerSize', 5)
hold off
xlabel 'N_{Na}'
ylabel 'RS'
legend('SSA', 'RSDE', 'SDE')

% Fig. 10 in Dangerfield et al. Phys Rev E 2012

%Plotting Differences in RS
diff_den_r=abs(RS_ssa_d-RS_rsde_d);
diff_a_r=abs(RS_ssa_a-RS_rsde_a);
diff_den_s=abs(RS_ssa_d-RS_sde_d);
diff_a_s=abs(RS_ssa_a-RS_sde_a);

% Stimulus amplitude - values taken from findIthRS.m available at http://www.ece.mcmaster.ca/~ibruce/
% Stimulus amplitudes for constant area simulations

for i=1:length(I_amp(3, :))
    if isnan(jitter_ssa_a(3, i))
        jit_ssa(i)=0;
    else
        jit_ssa(i)=jitter_ssa_a(3, i);
    end
    if isnan(jitter_rsde_a(3, i))
        jit_rsde(i)=0;
    else
        jit_rsde(i)=jitter_rsde_a(3, i);
    end
    if isnan(jitter_sde_a(3, i))
        jit_sde(i)=0;
    else
        jit_sde(i)=jitter_sde_a(3, i);
    end
    if isnan(meanlt_ssa_a(3, i))
        ml_ssa(i)=0;
    else
        ml_ssa(i)=meanlt_ssa_a(3, i);
    end
    if isnan(meanlt_rsde_a(3, i))
        ml_rsde(i)=0;
    else
        ml_rsde(i)=meanlt_rsde_a(3, i);
    end
    if isnan(meanlt_sde_a(3, i))
        ml_sde(i)=0;
    else
        ml_sde(i)=meanlt_sde_a(3, i);
    end
    
end

diff_jitter_r=abs(jit_ssa-jit_rsde);
diff_jitter_s=abs(jit_ssa-jit_sde);
diff_ml_r=abs(ml_ssa-ml_rsde);
diff_ml_s=abs(ml_ssa-ml_sde);

figure
%title ('constant density')
subplot(2, 2, 1)
semilogx(Nna, diff_den_r, 'r', 'LineWidth', 3,'MarkerSize', 5)
hold on
semilogx(Nna, diff_den_s, 'o--b', 'LineWidth', 3,'MarkerSize', 5)
hold off
xlabel 'N_{Na}'
ylabel 'Difference in RS'
%title ('constant area')
subplot(2, 2, 2)
semilogx(Nna, diff_a_r, 'r', 'LineWidth', 3,'MarkerSize', 5)
hold on
semilogx(Nna, diff_a_s, 'o--b', 'LineWidth', 3,'MarkerSize', 5)
hold off
xlabel 'N_{Na}'
ylabel 'Difference in RS'
subplot(2, 2, 3)
plot(I_amp(3, :), diff_jitter_r, 'r', 'LineWidth', 3,'MarkerSize', 5)
hold on
plot(I_amp(3, :), diff_jitter_s, 'o--b', 'LineWidth', 3,'MarkerSize', 5)
hold off
xlabel 'I'
ylabel 'Difference in Jitter'
subplot(2, 2, 4)
plot(I_amp(3, :), diff_ml_r, 'r', 'LineWidth', 3,'MarkerSize', 5)
hold on
plot(I_amp(3, :), diff_ml_s, 'o--b', 'LineWidth', 3,'MarkerSize', 5)
hold off
xlabel 'I'
ylabel 'Difference in Latency'
legend ('RSDE', 'SDE')



function [xprime, t]=mh_gate_HH_ODE(t, x, V)

alpham=(1.872*(V-25.41))/(1-exp((25.410-V)/6.06)); %from Bruce 2009 paper
betam=(3.973*(21.001-V))/(1-exp((V-21.001)/9.41)); %from Bruce 2009 paper
alphah=(-0.549*(27.74+V))/(1-exp((V+27.74)/9.06)); %from Bruce 2009 paper
betah=22.57/(1+exp((56.0-V)/12.5)); %from Bruce 2009 paper

xprime=zeros(8, 1);

M=[-(betah+3*betam) alpham 0 0 alphah 0 0 0; 3*betam -(alpham+2*betam+betah) 2*alpham 0 0 alphah 0 0; 0 2*betam -(2*alpham+betam+betah) 3*alpham 0 0 alphah 0; 0 0 betam -(3*alpham+betah) 0 0 0 alphah];
M=[M; betah 0 0 0 -(alphah+3*betam) alpham 0 0; 0 betah 0 0 3*betam -(alphah+alpham+2*betam) 2*alpham 0; 0 0 betah 0 0 2*betam -(2*alpham+betam+alphah) 3*alpham; 0 0 0 betah 0 0 betam -(alphah+3*alpham)];

% ODE for Markov model of sodium channel
xprime=M*x;

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
    
    X(1, i)=X(1, i-1)+(1/C*((-X(1, i-1)/R)-g*round(N*X(2, i-1))*(X(1, i-1)-Ena)+Istim(i-1)))*Dt;
    
    % Update
    X(2, i)=na_channel(1, i);
    t(i)=t(i-1)+Dt;
    
end

function [X, t] = HH_SDE_FE_sims (t_fin, X0, na0, Dt, I_amp, N, scaled)

%%% Code adapted from code for (Goldwyn et al., PLoS Comp Bio 2011)
%%% available on ModelDB website, accession number 138950

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
    X(1, i)=X(1, i-1)+(1/C*((-X(1, i-1)/R)-g*round(N*X(2, i-1))*(X(1, i-1)-Ena)+Istim(i-1)))*Dt;
    
    % Update
    X(2, i)=na_channel(1, i);
    t(i)=t(i-1)+Dt;
    
end

function err=intgausserr(X0,x,y)
% err=intgausserr(X0,x,y)
%
% Provides error to fit of Eq. (16) in:
%
%    Bruce, I. C. (2009). "Evaluation of stochastic differential equation
%    approximation of ion channel gating models," Annals of Biomedical
%    Engineering 37(4):824�838.

%%% This code is taken directly from the authors website and is available at
%%% http://www.ece.mcmaster.ca/~ibruce/

thrsh    = X0(1);
noisestd = X0(2);

f = 1/2*(erf((x-thrsh)/sqrt(2)./noisestd)+1);

err = norm(f-y);
% Ian C. Bruce (ibruce@ieee.org) � 2009

function x = projsplx(y)
% project an n-dim vector y to the simplex Dn
% Dn = { x : x n-dim, 1 >= x >= 0, sum(x) = 1}

% (c) Xiaojing Ye
% xyex19@gmail.com
%
% Algorithm is explained as in the linked document
% http://arxiv.org/abs/1101.6081
% or
% http://ufdc.ufl.edu/IR00000353/
%
% Jan. 14, 2011.

m = length(y); bget = false;

s = sort(y,'descend'); tmpsum = 0;

for ii = 1:m-1
    tmpsum = tmpsum + s(ii);
    tmax = (tmpsum - 1)/ii;
    if tmax >= s(ii+1)
        bget = true;
        break;
    end
end
    
if ~bget, tmax = (tmpsum + s(m) -1)/m; end;

x = max(y-tmax,0);

return;