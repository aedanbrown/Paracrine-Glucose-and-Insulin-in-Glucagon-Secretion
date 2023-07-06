
% This file runs the complete alpha- and beta-cell model in a perfusion
% setting with different step changes in glucose in flow rates in an "in
% vivo" setting. The perfusion rate, volume, and
% number of islets are set to in vivo values. The glucose in concentration
% is set to 5 mM. Instead of a step change in insulin and glucagon in flow
% concentration, it is fixed, which better reflects the in vivo setting.

%Experimental conditions
numIslets = 10^6; %Number of islets
g_in_0 = 5; %mM - Initial glucose in flow concentration
betaOnOff = true; %true = on, false = off
alphaOnOff = true; %true = on, false = off
I_t_in_c = 2.2e-5; %mg/dL - Insulin in flow step concentration
G_t_in_c = 2.5e-6; %mg/dL - Glucagon in flow step concentration

g_up_vals =  5; %mM - range of glucose values to examine


%Parameters
%--------------------------------------------------------------------------
%Basal levels for humans
gba_ = 90.08; %mg/dL - basal glucose levels (Alcazar & Buchwald, 2019)
Gba_ = 2.5e-6; %mg/dL - basal glucagon levels (van Vliet et. al., 2020 and http://www.ncbi.nlm.nih.gov/books/NBK279127)
Iba_ = 2.2e-5; %mg/dL - basal insulin levels (Van Vliet et. al., 2020)

%--------------------------------------------------------------------------
%Kinetic secretion parameters
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
%Steady-state model parameters

%Interaction parameters from mice

%Net beta-cell signal
m_GB_ = 1.11; %[] - maximum contribution of glucagon to beta cell signal
h_GB_ = 502; %[] - half-maximum value for contribution of glucagon
n_GB_ = 0.628; %[] - Hill exponent for contribution of insulin
h_gB_ = 1.07; %[] - half-maximum value for on/off effect of glucose
n_gB_ = 0.350; %[] - Hill exponent for on/off effect of glucose

%Net alpha-cell signal
h_IA_ = 10.0; %[] - half-maximum for contribution of insulin to alpha cell signal
n_IA_ = 1.17; %[] - Hill exponent for contribution of insulin 
m_g_ = 0.600; %[] Maximum fraction of glucose signal that insulin can remove

%--------------------------------------------------------------------------

%Mass balance related parameters
Q_ = 0.0117; %dL/min - perfusion rate
V_P_ = 1; %dL - volume

%--------------------------------------------------------------------------
%Steady-state interaction parameters

%Steady-state insulin secretion
m_I_ = 1.03e-7./15.*numIslets.*betaOnOff; %mg/min - maximum steady-state insulin secretion rate
h_I_ = 3.97; %[] - half-maximum for insulin secretion rate
n_I_ = 4.84; %[] - Hill exponent for insulin secretion rate

%Steady-state glucagon secretion
m_G_ = 2.24e-9./15.*numIslets.*alphaOnOff; %mg/min - maximum steady-state glucagon secretion rate
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
      
%Step functions are not used, because this doesn't reflect in vivo
G_t_in = @(t) G_t_in_c.*ones(length(t),1); %Glucagon in flow concentration step function
I_t_in = @(t) I_t_in_c.*ones(length(t),1); %Insulin in flow concentration step function
t = 0:0.1:60; %min - time range to examine

n = length(g_up_vals);

R_I = zeros(length(g_up_vals),1);
R_G = zeros(length(g_up_vals),1);

for i = 1:n

    g_t_in = @(t) g_in_0.*18.016 + (t > 0).*(g_up_vals(i) - g_in_0).*18.016;
    [~,results] = simulate_alphaBetaModel_perfusion(params_,g_t_in,G_t_in,I_t_in,t);
    
    R_I(i) = (results(end,1) - results(1,1)).*V_P_ - Q_.*trapz(t,I_t_in(t)) + Q_.*trapz(t,results(:,1));
        %Because I_t_in is a function of time, need to calculate insulin
        %secretion by integrating the mass balance
    R_G(i) = (results(end,2) - results(1,2)).*V_P_ - Q_.*trapz(t,G_t_in(t)) + Q_.*trapz(t,results(:,2));
        %Because G_t_in is a function of time, need to calculate glucagon
        %secretion by integrating the mass balance
    
end


fprintf("Done\n")