function [mn_na_f, var_na_f, mn_k_f, var_k_f, mn_na_v, var_na_v, mn_k_v, var_k_v] = mean_variance_plots ();

%%% OUTPUTS
% mn_na_f = matrix of means for RSDE and SSA for sodium channel for range of
% channel number and range of voltage values (fixed in simulation)
%
% var_na_f = matrix of variances for RSDE and SSA for sodium channel for range of
% channel number and range of voltage values (fixed in simulation)
%
% mn_k_f = matrix of means for RSDE and SSA for potassium channel for range of
% channel number and range of voltage values (fixed in simulation)
%
% var_k_f = matrix of variances for RSDE and SSA for potassium channel for range of
% channel number and range of voltage values (fixed in simulation)
%
% mn_na_v = matrix of means for RSDE and SSA for sodium channel for range of
% channel number and varying voltage (fixed voltage path)
%
% var_na_v = matrix of variances for RSDE and SSA for sodium channel for range of
% channel number and varying voltage (fixed voltage path)
%
% mn_k_v = matrix of means for RSDE and SSA for potassium channel for range of
% channel number and varying voltage (fixed voltage path)
%
% var_k_v = matrix of variances for RSDE and SSA for potassium channel for range of
% channel number and varying voltage (fixed voltage path)
%
% Figures 4, 5, 6, 7 and 8 in Dangerfield et al. Phys Rev E 2012, Vol. 85 Issue 5 2012

clear all
close all
clc

% Number of simulations
Nsim=10;

% Total number of channels
Nna=[500 750 1000];
Nk=[500 750 1000];

%% FIXED VOLTAGE

V=[-45 -25 0 20 40 60]; % Voltage values
t_fin=100; % Length of simulation
Dt=0.01; % Time step used in simulation

% Initial conditions for RSDE
X0na=1/8*ones(8,1); % Sodium channel
X0k=1/5*ones(5,1); % Potassium channel

% Mean and variance matrices
% Top row is for RSDE middle is for SSA (for each channel number)
mn_na_f=zeros(2*length(Nna), length(V));
var_na_f=zeros(2*length(Nna), length(V));
mn_k_f=zeros(2*length(Nna), length(V));
var_k_f=zeros(2*length(Nna), length(V));

%%% MAIN LOOP
for i=1:length(Nna)
    
    % Set up initial condition for SSA simulations
    
    poss_na=round(Nna(i)*X0na(1:7));
    poss_na=[poss_na; Nna(i)-sum(poss_na)];
    % Need to ensure that sum of states is equal to total number of channels
    if min(poss_na)>0
        X0na_ssa=poss_na;
    else
        X0na_ssa=[round(Nna(i)/2); round(Nna(i)/2); Nna(i)-(round(Nna(i)/2)+round(Nna(i)/2)); 0; 0; 0; 0; 0];
    end
    
    poss_k=round(Nk(i)*X0k(1:4));
    poss_k=[poss_k; Nk(i)-sum(poss_k)];
    % Need to ensure that sum of states is equal to total number of channels
    if min(poss_k)>0
        X0k_ssa=poss_k;
    else
        X0k_ssa=[round(Nk(i)/2); round(Nk(i)/2); Nk(i)-(round(Nk(i)/2)+round(Nk(i)/2)); 0; 0];
    end
    
    for l=1:length(V)
        
        % sum of proportion of open channels for RSDE method
        mp_na=0;
        mp_k=0;
        
        % sum of squares of proportion of open channels for RSDE method
        sqp_na=0;
        sqp_k=0;
        
        % sum of number of open channels for SSA method
        mssa_na=0;
        mssa_k=0;
        
        % sum of squares of number of opend channels for SSA method
        sqssa_na=0;
        sqssa_k=0;
        
        for k=1:Nsim
            i
            l
            k
            
            % RSDE METHOD
            [X_na, t_na] = RSDE_Na_channel_fixed_v(t_fin, X0na, Dt, V(l), Nna(i));
            [X_k, t_k] = RSDE_Nk_channel_fixed_v(t_fin, X0k, Dt, V(l), Nk(i));
            % update sums
            mp_na=mp_na+X_na(1, end);
            mp_k=mp_k+X_k(1, end);
            sqp_na=sqp_na+X_na(1, end).^2;
            sqp_k=sqp_k+X_k(1, end).^2;
            
            % SSA METHOD
            [X_na_ssa,t_na_ssa]=SSA_Na_channel_fixed_v(t_fin, X0na_ssa, Dt, V(l));
            [X_k_ssa,t_k_ssa]=SSA_Nk_channel_fixed_v(t_fin, X0k_ssa, Dt, V(l));
            
            % update sums
            mssa_na=mssa_na+X_na_ssa(1, end);
            mssa_k=mssa_k+X_k_ssa(1, end);
            sqssa_na=sqssa_na+X_na_ssa(1, end).^2;
            sqssa_k=sqssa_k+X_k_ssa(1, end).^2;
            
        end
        frm=2*i-1;
        to=2*i;
        
        % Update mean and variance matrices
        mn_na_f(frm:to, l)=[mp_na; mssa_na/Nna(i)]./Nsim;
        mn_k_f(frm:to, l)=[mp_k; mssa_k/Nk(i)]./Nsim;
        var_na_f(frm:to, l)=[sqp_na; sqssa_na/(Nna(i)^2)]./Nsim-mn_na_f(frm:to, l).^2;
        var_k_f(frm:to, l)=[sqp_k; sqssa_k/(Nk(i)^2)]./Nsim-mn_k_f(frm:to, l).^2;
        
    end
    
end

%%% FIGURES 

% Mean and variance plots for N=1000
frm=5;
too=6;

% Figure 4 in Dangerfield et al. Phys Rev E 2012
figure('Name', 'Mean and Variance Plots as a function of V for 1000 channels')
subplot(2, 2, 1)
plot(V, mn_na_f(frm, :), '-ob', 'LineWidth', 3,'MarkerSize', 5)
hold on
plot(V, mn_na_f(too, :), '+r','MarkerSize', 5)
hold off
xlabel 'Voltage'
ylabel 'Mean'
xlim([min(V), max(V)])
subplot(2, 2, 2)
plot(V, sqrt(var_na_f(frm, :)), '-ob', 'LineWidth', 3,'MarkerSize', 5)
hold on
plot(V, sqrt(var_na_f(too, :)),'+r','MarkerSize', 5)
hold off
xlabel 'Voltage'
ylabel 'Standard Deviation'
xlim([min(V), max(V)])
subplot(2, 2, 3)
plot(V, mn_k_f(frm, :), '-ob','LineWidth', 3,'MarkerSize', 5)
hold on
plot(V, mn_k_f(too, :),'+r','MarkerSize', 5)
hold off
xlabel 'Voltage'
ylabel 'Mean'
xlim([min(V), max(V)])
subplot(2, 2, 4)
plot(V, sqrt(var_k_f(frm, :)), '-ob','LineWidth', 3,'MarkerSize', 5)
hold on
plot(V, sqrt(var_k_f(too, :)),'+r','MarkerSize', 5)
hold off
xlim([min(V), max(V)])
xlabel 'Voltage'
ylabel 'Standard Deviation'
legend('RSDE', 'SSA')

% Difference in mean and variance plots

% Differences between the mean and variance of SSA and RSDE
di_mn_na=zeros(length(Nna), length(V));
di_std_na=zeros(length(Nna), length(V));
di_mn_k=zeros(length(Nk), length(V));
di_std_k=zeros(length(Nk), length(V));
for i=1:length(Nna)
    frm=2*i-1;
    too=2*i;
    di_mn_na(i, :)=abs(mn_na_f(frm,:)-mn_na_f(too, :));
    di_std_na(i, :)=abs(sqrt(var_na_f(frm,:))-sqrt(var_na_f(too, :)));
end
for i=1:length(Nk)
    frm=2*i-1;
    too=2*i;
    di_mn_k(i, :)=abs(mn_k_f(frm,:)-mn_k_f(too, :));
    di_std_k(i, :)=abs(sqrt(var_k_f(frm,:))-sqrt(var_k_f(too, :)));
end

% Figure 5 in Dangerfield et al. Phys Rev E 2012
figure('Name', 'Difference Plots for Fixed Voltage')
subplot(2, 2, 1)
plot(V, di_mn_na(1, :), 'b', 'LineWidth', 3)
hold on
plot(V, di_mn_na(2, :), '.k', 'LineWidth', 3, 'MarkerSize', 6)
plot(V, di_mn_na(3, :), '--r', 'LineWidth', 3)
hold off
xlabel 'Voltage'
ylabel 'Difference in Means'
xlim([min(V), max(V)])
subplot(2, 2, 2)
plot(V, di_std_na(1, :), 'b', 'LineWidth', 3)
hold on
plot(V, di_std_na(2, :), '.k', 'LineWidth', 3, 'MarkerSize', 6)
plot(V, di_std_na(3, :), '--r', 'LineWidth', 3)
hold off
xlabel 'Voltage'
ylabel 'Difference in Std'
xlim([min(V), max(V)])
subplot(2, 2, 3)
plot(V, di_mn_k(1, :), 'b', 'LineWidth', 3)
hold on
plot(V, di_mn_k(2, :), '.k', 'LineWidth', 3, 'MarkerSize', 6)
plot(V, di_mn_k(3, :), '--r', 'LineWidth', 3)
hold off
xlabel 'Voltage'
ylabel 'Difference in Means'
xlim([min(V), max(V)])
subplot(2, 2, 4)
plot(V, di_std_k(1, :), 'b', 'LineWidth', 3)
hold on
plot(V, di_std_k(2, :), '.k', 'LineWidth', 3, 'MarkerSize', 6)
plot(V, di_std_k(3, :), '--r', 'LineWidth', 3)
xlabel 'Voltage'
ylabel 'Difference in Std'
xlim([min(V), max(V)])
legend ('N=500', 'N=750', 'N=1000')

% Clear variables used for varying voltage simulations

clear t_fin Dt V X0k X0na X0na_ssa X0k_ssa mp_na mp_k mssa_na mssa_k sqp_na sqp_k sqssa_na sqssa_k di_mn_na di_mn_k di_std_na di_std_k

%% VARYING VOLTAGE (FIXED VOLTAGE PATH)

t_fin=10; % Length of simulation
Dt=0.01; % Time step used in simulation
Nstep = ceil(t_fin/Dt);

% Set up voltage path by solving deterministic Dodgkin-Huxley model over [0 t_fin]
X0det=[0 0.5 0.5 0.5]';
tspan=0:Dt:t_fin;
[t, y]=ode45(@hodgkin_huxley_orig,tspan,X0det);
Vpath=y(:, 1); % Fixed voltage path used in stochastic simulation

% Set up initial condition - start system at steady state
V=0;
alpham=0.1*((25-V)/(exp((25-V)/10)-1));
betam=4*exp(-V/18);
alphah=0.07*exp(-V/20);
betah=1/(exp((30-V)/10)+1);
alphan=(0.01*(10-V))/(exp((10-V)/10)-1);
betan=0.125*exp(-V/80);

X0k=(1/(alphan+betan)^4)*[alphan^4*nchoosek(4, 4)
    alphan^3*betan*nchoosek(4, 3)
    alphan^2*betan^2*nchoosek(4, 2)
    alphan*betan^3*nchoosek(4, 1)
    betan^4*nchoosek(4, 0)];
X0na=1/((alpham+betam)^3*(alphah+betah))*[alpham^3*alphah*nchoosek(3,3)
    alpham^2*betam*alphah*nchoosek(3,2)
    alpham*betam^2*alphah*nchoosek(3,1)
    betam^3*alphah*nchoosek(3,0)
    alpham^3*betah*nchoosek(3,3)
    alpham^2*betam*betah*nchoosek(3,2)
    alpham*betam^2*betah*nchoosek(3,1)
    betam^3*betah*nchoosek(3,0)];

% Mean and variance matrices
% Top row is for RSDE middle is for SSA (for each channel number)
mn_na_v=zeros(2*length(Nna), Nstep+1);
mn_k_v=zeros(2*length(Nk), Nstep+1);
var_na_v=zeros(2*length(Nna), Nstep+1);
var_k_v=zeros(2*length(Nk), Nstep+1);

% MAIN LOOP
for i=1:length(Nna)
    i
    
    % Set initial condition for SSA 
    x0_poss=round(Nna(i)*X0na(1:7));
    x0_poss2=[x0_poss ; Nna(i)-sum(x0_poss)];
    % Need to ensure that sum of states is equal to total number of channels
    if min(x0_poss2)>=0
        X0na_ssa=x0_poss2;
    else
        ex=sum(x0_poss)-Nna(i);
        Ind=find(x0_poss>ex);
        if Ind(1)==1
            X0na_ssa=[x0_poss(1)-ex; x0_poss(2:7) ; 0];
        elseif Ind(1)==7
            X0na_ssa=[x0_poss(1:Ind(1)-1); x0_poss(Ind(1))-ex; 0];
        else
            X0na_ssa=[x0_poss(1:Ind(1)-1); x0_poss(Ind(1))-ex; x0_poss(Ind(1)+1:end); 0];
        end
    end
    poss_na=round(Nna(i)*X0na(1:7));
    poss_na=[poss_na; Nna(i)-sum(poss_na)];
    % Need to ensure that sum of states is equal to total number of channels
    if min(poss_na)>=0
        X0na_ssa=poss_na;
    else
        error('SSA condition sum too large for Na')
    end
    y0_poss=round(Nk(i)*X0k(1:4));
    y0_poss2=[y0_poss ; Nk(i)-sum(y0_poss)];
    if min(y0_poss2)>=0
        X0k_ssa=y0_poss2;
    else
        ex=sum(y0_poss)-Nk(i);
        Ind=find(y0_poss>ex);
        if Ind(1)==1
            X0k_ssa=[y0_poss(1)-ex; y0_poss(2:4) ; 0];
        elseif Ind(1)==4
            X0k_ssa=[y0_poss(1:Ind(1)-1); y0_poss(Ind(1))-ex; 0];
        else
            X0k_ssa=[y0_poss(1:Ind(1)-1); y0_poss(Ind(1))-ex; y0_poss(Ind(1)+1:end); 0];
        end
    end
    poss_k=round(Nk(i)*X0k(1:4));
    poss_k=[poss_k; Nk(i)-sum(poss_k)];
    if min(poss_k)>0
        X0k_ssa=poss_k;
    else
         error('SSA condition sum too large for k')
    end
    
    
    % sum of proportion of open channels for RSDE method
    mp_na=zeros(1, Nstep+1);
    mp_k=zeros(1, Nstep+1);
    
    % sum of squares of proportion of open channels for RSDE method
    sqp_na=zeros(1, Nstep+1);
    sqp_k=zeros(1, Nstep+1);
    
    % sum of number of open channels for SSA method
    mssa_na=zeros(1, Nstep+1);
    mssa_k=zeros(1, Nstep+1);
    
     % sum of squares of number of open channels for SSA method
    sqssa_na=zeros(1, Nstep+1);
    sqssa_k=zeros(1, Nstep+1);
        
    for k=1:Nsim
         i
         
         k
        % RSDE METHOD 
        [X_na, t_na] = RSDE_Na_channel_varying_v(t_fin, X0na, Dt, Nna(i), Vpath);
        [X_k, t_k] = RSDE_Nk_channel_varying_v(t_fin, X0k, Dt, Nk(i), Vpath);
        % Update sums
        mp_na=mp_na+X_na(1, :);
        mp_k=mp_k+X_k(1, :);
        sqp_na=sqp_na+X_na(1, :).^2;
        sqp_k=sqp_k+X_k(1, :).^2;
        
        % SSA METHOD
        [X_na_ssa, t_na_ssa] = SSA_Na_channel_varying_v(t_fin, X0na_ssa, Dt, Vpath);
        [X_k_ssa, t_k_ssa] = SSA_Nk_channel_varying_v(t_fin, X0k_ssa, Dt, Vpath);
        % Update sums
        mssa_na=mssa_na+X_na_ssa(1, :);
        mssa_k=mssa_k+X_k_ssa(1, :);
        sqssa_na=sqssa_na+X_na_ssa(1, :).^2;
        sqssa_k=sqssa_k+X_k_ssa(1, :).^2;
        
    end
    
    frm=2*i-1;
    to=2*i;
    
    % Update mean and variance matrices
    mn_na_v(frm:to, :)=[mp_na; mssa_na/Nna(i)]./Nsim;
    mn_k_v(frm:to, :)=[mp_k; mssa_k/Nk(i)]./Nsim;
    var_na_v(frm:to, :)=[sqp_na; sqssa_na/Nna(i)^2]./Nsim-mn_na_v(frm:to, :).^2;
    var_k_v(frm:to, :)=[sqp_k; sqssa_k/Nk(i)^2]./Nsim-mn_k_v(frm:to, :).^2;
  
end

%%% FIGURES

% Figure 6 in Dangerfield et al. Phys Rev E 2012
% Plot of the voltage path used
figure('Name', 'Voltage Path')
plot(t, Vpath)
xlabel 'Time'
ylabel 'Voltage'

% Mean and variance plots for number of channels equal to 1000

frm=5;
too=6;

% Figure 7 in Dangerfield et al. Phys Rev E 2012
figure('Name', 'Mean and Variance Plots as a function of time for 1000 channels')
subplot(2, 2, 1)
plot(t, mn_na_v(frm, :), '-ob', 'LineWidth', 3,'MarkerSize', 5)
hold on
plot(t, mn_na_v(too, :), '+r','MarkerSize', 5)
hold off
xlabel 'Time'
ylabel 'Mean'
subplot(2, 2, 2)
plot(t(2:end), sqrt(var_na_v(frm, 2:end)), '-ob', 'LineWidth', 3,'MarkerSize', 5)
hold on
plot(t(2:end), sqrt(var_na_v(too, 2:end)), '+r','MarkerSize', 5)
hold off
xlabel 'Time'
ylabel 'Std'
subplot(2, 2, 3)
plot(t, mn_k_v(frm, :), '-ob', 'LineWidth', 3,'MarkerSize', 5)
hold on
plot(t, mn_k_v(too, :), '+r','MarkerSize', 5)
hold off
xlabel 'Time'
ylabel 'Mean'
subplot(2, 2, 4)
plot(t(2:end), sqrt(var_k_v(frm, 2:end)), '-ob', 'LineWidth', 3,'MarkerSize', 5)
hold on
plot(t(2:end), sqrt(var_k_v(too, 2:end)), '+r','MarkerSize', 5)
hold off
xlabel 'Time'
ylabel 'Std'
legend('RSDE', 'SSA')
    

% Differences in mean and variance between RSDE and SSA plots

%plotting differences in the means and variances
di_mn_na=zeros(length(Nna), length(t));
di_std_na=zeros(length(Nna), length(t));
di_mn_k=zeros(length(Nna), length(t));
di_std_k=zeros(length(Nna), length(t));

for i=1:length(Nna)
    frm=2*i-1;
    too=2*i;
    di_mn_na(i, :)=abs((mn_na_v(frm,:)-mn_na_v(too, :)));
    di_std_na(i, :)=abs((sqrt(var_na_v(frm,:))-sqrt(var_na_v(too, :))));
    di_mn_k(i, :)=abs((mn_k_v(frm,:)-mn_k_v(too, :)));
    di_std_k(i, :)=abs((sqrt(var_k_v(frm,:))-sqrt(var_k_v(too, :))));
end

% Figure 8 in Dangerfield et al. Phys Rev E 2012
figure('Name', 'Difference Plots for Varying Voltage')
subplot(2, 2, 1)
plot(t, di_mn_na(1, :), 'b', 'LineWidth', 3)
hold on
plot(t, di_mn_na(2, :), '.k', 'LineWidth', 3, 'MarkerSize', 6)
plot(t, di_mn_na(3, :), '--r', 'LineWidth', 3)
hold off
xlabel 'Time'
ylabel 'Difference in Means'
xlim([min(t), max(t)])
subplot(2, 2, 2)
plot(t, di_std_na(1, :), 'b', 'LineWidth', 3)
hold on
plot(t, di_std_na(2, :), '.k','LineWidth', 3, 'MarkerSize', 6)
plot(t, di_std_na(3, :), '--r', 'LineWidth', 3)
hold off
xlabel 'Time'
ylabel 'Difference in StD'
xlim([min(t), max(t)])
subplot(2, 2, 3)
plot(t, di_mn_k(1, :), 'b', 'LineWidth', 3)
hold on
plot(t, di_mn_k(2, :), '.k', 'LineWidth', 3, 'MarkerSize', 6)
plot(t, di_mn_k(3, :), '--r', 'LineWidth', 3)
hold off
xlabel 'Voltage'
ylabel 'Difference in Means'
xlim([min(t), max(t)])
subplot(2, 2, 4)
plot(t, di_std_k(1, :), 'b', 'LineWidth', 3)
hold on
plot(t, di_std_k(2, :), '.k', 'LineWidth', 3, 'MarkerSize', 6)
plot(t, di_std_k(3, :), '--r', 'LineWidth', 3)
hold off
xlabel 'Voltage'
ylabel 'Difference in StD'
xlim([min(t), max(t)])
legend ('N=500', 'N=750', 'N=1000')

%% FUNCTIONS USED IN SIMULATIONS 

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
function [X,t]=SSA_Nk_channel_fixed_v(t_fin, X0, Dt, V)

%%% Code adapted from code for (Bruce, ABME 2009) available at
%%% http://www.ece.mcmaster.ca/~ibruce/ (authors website)

%%% INPUTS
% t_fin is final time of solution (ms)
% X0 is vector with NUMBER (so entries of X0 are integers) of channels in each state at beginning of
% the simulation (initial state of the channel), ordering of vector is as
% follows - X0=[n4 n3 n2 n1 n0]'
% Dt is the time step
% v is the voltage value, which remains fixed throughout the simulation

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

% Transition rates, as function of the voltage
alphan=(0.01*(10-V))/(exp((10-V)/10)-1);%original HH model
betan=0.125*exp(-V/80);%original HH model

zeta=[4*betan
        alphan+3*betan
        2*alphan+2*betan
        3*alphan+betan
        4*alphan]; % Eq. (A3) of Mino et al. (ABME 2002)

% Main loop
for i=2:Nstep+1
  
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

function [X, t] = RSDE_Na_channel_varying_v(t_fin, X0, Dt, N, Vpath)

%%% INPUTS
% t_fin is final time of solution (ms)
% X0 is vector with proportion of channels in each state at beginning of
% the simulation (initial state of the channel), ordering of vector is as
% follows - X0=[m3h1 m2h1 m1h1 m0h1 m3h0 m2h0 m1h0 m0h0]'
% Dt is the time step
% N is the number of channels
% Vpath is the voltage path (so V varies over time)

%%% OUTPUTS
% t is time
% X - X(:, i) is the proportion of channels in each state at time t(i) 
%(ordering of X(:, i) the salphame as X0).

% Initialise size of outputs
Nstep = ceil(t_fin/Dt);     % Number of steps of size Dt=T/N
X=zeros(8,Nstep+1); 
t=zeros(Nstep+1,1);

% Set up initial conditions
X(:, 1)=X0;
t(1)=0;

% Initialise matrices for the RSDE
E=[-1 0 0 0 0 0 0 0 0 -1
    1 -1 0 0 0 0 0 0 -1 0
    0 1 -1 0 0 0 0 -1 0 0
    0 0 1 0 0 0 -1 0 0 0
    0 0 0 -1 0 0 0 0 0 1
    0 0 0 1 -1 0 0 0 1 0
    0 0 0 0 1 -1 0 1 0 0
    0 0 0 0 0 1 1 0 0 0];

% Brownian increments (scaled with the time step)  
dW = sqrt(Dt)*randn(10,Nstep);    

% Main loop
for i=2:Nstep+1
    
    v=Vpath(i);
    
    % Calculate transition rates at current timestep, as function of the voltage
    alpham=0.1*((25-v)/(exp((25-v)/10)-1));
    betam=4*exp(-v/18);
    alphah=0.07*exp(-v/20);
    betah=1/(exp((30-v)/10)+1);
    
    % Calculate drift term
    A=[-(betah+3*betam) alpham 0 0 alphah 0 0 0; 3*betam -(alpham+2*betam+betah) 2*alpham 0 0 alphah 0 0; 0 2*betam -(2*alpham+betam+betah) 3*alpham 0 0 alphah 0; 0 0 betam -(3*alpham+betah) 0 0 0 alphah];
    A=[A; betah 0 0 0 -(alphah+3*betam) alpham 0 0; 0 betah 0 0 3*betam -(alphah+alpham+2*betam) 2*alpham 0; 0 0 betah 0 0 2*betam -(2*alpham+betam+alphah) 3*alpham; 0 0 0 betah 0 0 betam -(alphah+3*alpham)];
    
    % Calculate diffusion term
    S1=sqrt((alpham*X(2, i-1)+3*betam*X(1, i-1)));
    S2=sqrt((2*alpham*X(3, i-1)+2*betam*X(2, i-1)));
    S3=sqrt((3*alpham*X(4, i-1)+betam*X(3, i-1)));
    S4=sqrt((alpham*X(6, i-1)+3*betam*X(5, i-1)));
    S5=sqrt((2*alpham*X(7, i-1)+2*betam*X(6, i-1)));
    S6=sqrt((3*alpham*X(8, i-1)+betam*X(7, i-1)));
    S7=sqrt((betah*X(4, i-1)+alphah*X(8, i-1)));
    S8=sqrt((betah*X(3, i-1)+alphah*X(7, i-1)));
    S9=sqrt((betah*X(2, i-1)+alphah*X(6, i-1)));
    S10=sqrt((betah*X(1, i-1)+alphah*X(5, i-1)));
    Noi=[S1*dW(1, i-1); S2*dW(2, i-1); S3*dW(3, i-1); S4*dW(4, i-1); S5*dW(5, i-1); S6*dW(6, i-1); S7*dW(7, i-1); S8*dW(8, i-1); S9*dW(9, i-1); S10*dW(10, i-1)];
    
    % Calculate SDE without reflecting term, at next time step
    Xnext=X(:, i-1)+A*X(:, i-1)*Dt+(1/sqrt(N))*E*Noi;
    
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

function [X,t]=SSA_Na_channel_varying_v(t_fin, X0, Dt, Vpath)

%%% Code adapted from code for (Bruce, ABME 2009) available at
%%% http://www.ece.mcmaster.ca/~ibruce/ (authors website)

%%% INPUTS
% t_fin is final time of solution (ms)
% X0 is vector with NUMBER (so entries of X0 are integers) of channels in each state at beginning of
% the simulation (initial state of the channel), ordering of vector is as
% follows - X0=[m3h1 m2h1 m1h1 m0h1 m3h0 m2h0 m1h0 m0h0]'
% Dt is the time step
% Vpath is the voltage path (so V varies over time)

%%% OUTPUTS
% t is time
% X - X(:, i) is the NUMBER of channels in each state at time t(i) 
%(ordering of X(:, i) the same as X0).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise size of outputs
Nstep = ceil(t_fin/Dt); 
X=zeros(8, Nstep+1); %stores the voltage and proportion of open sodium channels at time intervals of length Dt.
t=zeros(Nstep+1,1);

% Vectors describing state transitions from current to next state
curr_state=[0 4 3 2 8 7 6 1 2 3 5 6 7 8 7 6 5 4 3 2 1];
next_state=[0 3 2 1 7 6 5 2 3 4 6 7 8 4 3 2 1 8 7 6 5];

% Set up initial conditions
X(:, 1)=X0;
mh=X(:, 1);
t(1)=0;

% Main loop
for i=2:Nstep+1
    
    v=Vpath(i-1); % voltage constnat over time step of length Dt
    
    % Transition rates, as function of the voltage
    alpham=0.1*((25-v)/(exp((25-v)/10)-1));
    betam=4*exp(-v/18);
    alphah=0.07*exp(-v/20);
    betah=1/(exp((30-v)/10)+1);
    
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


function yprime = hodgkin_huxley_orig(t, y)

%%% Hodgkin-Huxley model ODEs

% Parameter values
C=1;
GN=120;
GK=36;
GL=0.3;
EN=115;
EK=-12;
EL=10.6;
I=10;
yprime=zeros(4, 1);

v=y(1);
m=y(2);
h=y(3);
n=y(4);

alpham=0.1*((25-v)/(exp((25-v)/10)-1)); %orig HH
betam=4*exp(-v/18);%orig HH
alphah=0.07*exp(-v/20);%orig HH
betah=1/(exp((30-v)/10)+1);%orig HH
alphan=(0.01*(10-v))/(exp((10-v)/10)-1);%orig HH
betan=0.125*exp(-v/80);%orig HH

% Voltage equation
yprime(1)=(1/C*(-GN*m^3*h*(v-EN)-GK*n^4*(v-EK)-GL*(v-EL)+I));
% Gating variable equations
yprime(2)=alpham*(1-m)-betam*m;
yprime(3)=alphah*(1-h)-betah*h;
yprime(4)=alphan*(1-n)-betan*n;

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