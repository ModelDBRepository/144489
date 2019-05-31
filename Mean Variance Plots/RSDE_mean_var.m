%% Calculates the mean and variance of number of open channels for the

close all
clear all
clc

Nsim=100000;
Nna=[500 750 1000];
Nk=[500 750 1000];

%% FIXED VOLTAGE

V=[-45 -25 0 20 40 60];
t_fin=100;
Dt=0.01;

X0na=1/8*ones(8,1);
X0k=1/5*ones(5,1);

mn_na=zeros(length(Nna), length(V)); 
var_na=zeros(length(Nna), length(V));
mn_k=zeros(length(Nna), length(V)); 
var_k=zeros(length(Nna), length(V));

for i=1:length(Nna)
    
    for l=1:length(V)
        
        mp_na=0;
        mp_k=0;
        sqp_na=0;
        sqp_k=0;
        
        for k=1:Nsim
            i
            l
            k
            [X_na, t_na] = RSDE_Na_channel_fixed_v(t_fin, X0na, Dt, V(l), Nna(i));
            [X_k, t_k] = RSDE_Nk_channel_fixed_v(t_fin, X0k, Dt, V(l), Nk(i));
            mp_na=mp_na+X_na(1, end);
            mp_k=mp_k+X_k(1, end);
            sqp_na=sqp_na+X_na(1, end)^2;
            sqp_k=sqp_k+X_k(1, end)^2;
            
        end
        
        mn_na(i, l)=mp_na./Nsim;
        mn_k(i, l)=mp_k./Nsim;
        var_na(i, l)=sqp_na./Nsim-mn_na(i, l).^2;
        var_k(i, l)=sqp_k./Nsim-mn_k(i, l).^2;
        
    end
 
end
   figure
    subplot(2, 2, 1)
    plot(V, mn_na)
    subplot(2, 2, 2)
    plot(V, var_na)
    subplot(2, 2, 3)
    plot(V, mn_k)
    subplot(2, 2, 4)
    plot(V, var_k)
    legend ('500', '750', '10000')
%% VARYING VOLTAGE (FIXED VOLTAGE PATH)
t_fin=10;
Dt=0.01;
Nstep = ceil(t_fin/Dt);

X0det=[0 0.5 0.5 0.5]';
tspan=0:Dt:t_fin;
[t, y]=ode45(@hodgkin_huxley_orig,tspan,X0det);
Vpath=y(:, 1);

V=0;
alpham=0.1*((25-V)/(exp((25-V)/10)-1));%original HH model
betam=4*exp(-V/18);%original HH model
alphah=0.07*exp(-V/20);%original HH model
betah=1/(exp((30-V)/10)+1);%original HH model
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

mn_na=zeros(length(Nna), Nstep+1);
mn_nk=zeros(length(Nk), Nstep+1);
var_na=zeros(length(Nna), Nstep+1);
var_nk=zeros(length(Nk), Nstep+1);

for i=1:length(Nna)
    i
        
    mp_na=zeros(1, Nstep+1);
    mp_k=zeros(1, Nstep+1);
    sqp_na=zeros(1, Nstep+1);
    sqp_k=zeros(1, Nstep+1);
    
    
    for k=1:Nsim
         i
         l
         k
        [X_na, t_na] = RSDE_Na_channel_varying_v(t_fin, X0na, Dt, Nna(i), Vpath);
        [X_k, t_k] = RSDE_Nk_channel_varying_v(t_fin, X0k, Dt, Nk(i), Vpath);
        
        mp_na=mp_na+X_na(1, :);
        mp_k=mp_k+X_k(1, :);
        sqp_na=sqp_na+X_na(1, :).^2;
        sqp_k=sqp_k+X_k(1, :).^2;
        
    end
    
    mn_na(i, :)=mp_na./Nsim;
    mn_nk(i, :)=mp_k./Nsim;
    
    var_na(i, :)=sqp_na./Nsim-mn_na(i, :).^2;
    var_nk(i, :)=sqp_k./Nsim-mn_nk(i, :).^2;

    figure
    subplot(2, 2, 1)
    plot(t, mn_na(i, :))
    subplot(2, 2, 2)
    plot(t, var_na(i, :))
    subplot(2, 2, 3)
    plot(t, mn_nk(i, :))
    subplot(2, 2, 4)
    plot(t, var_nk(i, :))
    
end