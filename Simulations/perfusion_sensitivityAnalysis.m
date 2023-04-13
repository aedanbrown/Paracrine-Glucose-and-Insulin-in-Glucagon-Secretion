%This file runs the code used to perform sensitivity analysis.


%Experimental conditions
numIslets = 15;
g_in_0 = 1; %mM
betaOnOff = true; %true = on, false = off
alphaOnOff = true; %true = on, false = off
I_t_in_c = 0; %mg/dL
G_t_in_c = 0; %mg/dL

g_up_val = 15; %mM

%32 parameters * 2 measurements * 11 values = 706 array entries

%Parameters
%--------------------------------------------------------------------------
%Basal levels for humans
gba_ = 90.08; %mg/dL - basal glucose levels (Alcazar & Buchwald, 2019)
Gba_ = 2.5e-6; %mg/dL - basal glucagon levels (van Vliet et. al., 2020 and http://www.ncbi.nlm.nih.gov/books/NBK279127)
Iba_ = 2.2e-5; %mg/dL - basal insulin levels (Van Vliet et. al., 2020)

%--------------------------------------------------------------------------
%Kinetic secretion parameters - Curvefit to Henquin et. al. 2015 data
k_gB_ = 0.554; %1/min - beta cell glucose transduction
k_G_  = 0.554; %1/min - beta cell glucagon transduction
k_gA_ = 0.554/25; %1/min - alpha cell glucose transduction
k_I_  = 0.554*5; %1/min - alpha cell insulin transduction

%Rate constant for transfer from insulin pool 1 to insulin pool 2
m_I1_ = 0.336; %1/min - max value
h_I1_ = 3.75; %[] - half-maximum value
n_I1_ = 9.97; %[] - Hill exponent

%Rate constant for transfer from insulin pool 2 to outside cell
m_I2_ = 0.360; %1/min - max value
h_I2_ = 0.968; %[] - half-maximum value
n_I2_ = 6.68; %[] - Hill exponent

%Rate constant for transfer from glucagon pool 1 to glucagon pool 2
m_G1_ = 0.336; %1/min - max value
h_G1_ = 3.75; %[] - half-maximum value
n_G1_ = 9.97; %[] - Hill exponent

%Rate constant 
m_G2_ = 0.360; %1/min - max value
h_G2_ = 0.968; %[] - half-maximum value
n_G2_ = 6.68; %[] - Hill exponent


%--------------------------------------------------------------------------
%Steady-state model parameters - curvefit to data from Vieira et. al., 2007

%Interaction parameters from mice

%Net beta-cell signal
m_GB_ = 1.11; %[] - maximum contribution of glucagon to beta cell signal
h_GB_ = 502; %[] - half-maximum value for contribution of glucagon
n_GB_ = 0.628; %[] - Hill exponent for contribution of insulin
h_gB_ = 1.07; %[] - half-maximum value for on/off effect of glucose
n_gB_ = 0.350; %[] - Hill exponent for on/off effect of glucose

%Net alpha-cell signal
h_IA_ = 10.0; %[] - half-maximum for contribution of insulin to alpha cell signal - changed from mice value
n_IA_ = 1.17; %[] - Hill exponent for contribution of insulin - changed from mice value
m_g_ = 0.600; %[] Maximum fraction of glucose signal that insulin can remove

%--------------------------------------------------------------------------

%Mass balance related parameters
Q_ = 1/10^3*10; %dL/min - perfusion rate
V_P_ = 1/10^3*10; %dL - volume

%--------------------------------------------------------------------------
%Curvefit to Braun et. al., 2010 data - Braun et. al., 2010 used 10-20
%islets, so the secretion rate is divided by 15 islets

%Steady-state insulin secretion
m_I_ = 1.03e-7./15.*numIslets.*betaOnOff; %mg/min - maximum steady-state insulin secretion rate
h_I_ = 3.97; %[] - half-maximum for insulin secretion rate
n_I_ = 4.84; %[] - Hill exponent for insulin secretion rate

%Steady-state glucagon secretion
m_G_ = 2.24e-9./15.*numIslets.*alphaOnOff; %0.9e-8; %mg/min - maximum steady-state glucagon secretion rate
h_G_ = 1.06; %[] - half-maximum for glucagon secretion rate
n_G_ = 3.5; %[] - Hill exponent for glucagon secretion rate

%Background signals
X_B0_ = 2.60;
X_A0_ = 4.40;

%--------------------------------------------------------------------------


params_ = [gba_, Gba_, Iba_, ...
          k_gB_, k_G_, k_gA_, k_I_, ...
          m_GB_, h_GB_, n_GB_, h_gB_, n_gB_, X_B0_, ...
          h_IA_, n_IA_, X_A0_, m_g_, ...
          m_I_, h_I_, n_I_, ...
          m_G_, h_G_, n_G_, ...
          m_I1_, h_I1_, n_I1_, m_I2_, h_I2_, n_I2_, ...
          m_G1_, h_G1_, n_G1_, m_G2_, h_G2_, n_G2_, ...
          Q_,V_P_];
      
      
G_t_in = @(t) (t > 0).*G_t_in_c;
I_t_in = @(t) (t > 0).*I_t_in_c;
t = 0:0.1:60;


S_I_results = zeros(11,32);
S_G_results = zeros(11,32);
multipliers = [2/3 3/4 4/5 5/6 9/10 1 10/9 6/5 5/4 4/3 3/2];
    %Will multiply each parameter by this value
g_t_in = @(t) g_in_0.*18.016 + (t > 0).*(g_up_val - g_in_0).*18.016;

for i = 4:35 %First three parameters are basal levels
    
    params_adj = params_; %Reset the parameters
    
    for j = 1:11  
    
        params_adj(i) = params_(i).*multipliers(j); %Adjust the parameter of interest

        [~,y] = simulate_alphaBetaModel_perfusion(params_adj,g_t_in,G_t_in,I_t_in,t);
            %Run the simulation
    
        S_I_results(j,i-3) = (y(end,1) - y(1,1)).*V_P_ - Q_.*trapz(t,I_t_in(t)) + Q_.*trapz(t,y(:,1));
            %Because I_t_in is a function of time, need to calculate insulin
            %secretion by integrating the mass balance
        S_G_results(j,i-3) = (y(end,2) - y(1,2)).*V_P_ - Q_.*trapz(t,G_t_in(t)) + Q_.*trapz(t,y(:,2));
            %Because G_t_in is a function of time, need to calculate glucagon
            %secretion by integrating the mass balance

    end
    
end


fprintf("Done\n")