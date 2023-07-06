function dydt = alphaBetaModel_perfusion(t,y,p,g_t_in,G_t_in,I_t_in)

% alphabetaModel_perfusion runs the alpha- and beta-cell secretion model in
% a "perfusion" setting: there are in and out flow rates
% It takes in a scalar time t, a vector of the current state of the system 
% y, a vector of parameters p, a function handle g_t_in that represents
% the incoming glucose concentration as a function of time (can only take 
% in time), a function handle G_t_in that represents the incoming glucagon
% concentration as a function of time (can only take in time), and a
% function handle I_t_in that representes the incoming insulin
% concentration
% It returns the derivatives of insulin concentration, glucagon
% concentration, beta cell glucose signal, beta cell glucagon signal,
% insulin pool 1, insulin pool 2, alpha cell glucagon signal, alpha cell
% insulin signal, glucagon pool 1, glucagon pool 2, and glucose 
% concentration at the state described by the function inputs

    %Store parameters
    p = num2cell(p);
    [gba, Gba, Iba, ...
     k_gB, k_G, k_gA, k_I, ...
     m_GB, h_GB, n_GB, h_gB, n_gB, X_B0, ...
     h_IA, n_IA, X_A0, m_g, ...
     m_I, h_I, n_I, ...
     m_G, h_G, n_G, ...
     m_I1, h_I1, n_I1, m_I2, h_I2, n_I2, ...
     m_G1, h_G1, n_G1, m_G2, h_G2, n_G2, ...
     Q, V_P] = p{:};

    %Obtain current system values
    g = y(11); %Glucose concentration
    I = y(1); %Insulin concentration
    G = y(2); %Glucagon concentration

    X_gB = y(3); %Glucose signal in beta cells
    X_G = y(4); %Glucagon signal in beta cells
    I_1 = y(5); %Mass of insulin in first pool
    I_2 = y(6); %Mass of insulin in second pool

    X_gA = y(7); %Glucose signal in alpha cells
    X_I = y(8); %Insulin signal in alpha cells
    G_1 = y(9); %Mass of glucagon in first pool
    G_2 = y(10); %Mass of glucagon in second pool


    %Signal transduction
    dX_gB = k_gB.*(g./gba - X_gB);
    dX_G = k_G.*(G./Gba - X_G);
    dX_gA = k_gA.*(g./gba - X_gA);
    dX_I = k_I.*(I./Iba - X_I);

    %Net signals
    X_B = Y_B(X_gB,X_G,m_GB,h_GB,n_GB,h_gB,n_gB,X_B0);
    X_A = Y_A(X_gA,X_I,h_IA,n_IA,X_A0,m_g);

    %Steady-state secretion
    R_Iss_ = R_Iss(X_B,m_I,h_I,n_I);
    R_Gss_ = R_Gss(X_A,m_G,h_G,n_G);

    %Transient secretion
    dI_1 = R_Iss_ - hill(X_B,m_I1,h_I1,n_I1).*I_1;
    dI_2 = hill(X_B,m_I1,h_I1,n_I1).*I_1 - hill(X_B,m_I2,h_I2,n_I2).*I_2;
    R_I = hill(X_B,m_I2,h_I2,n_I2).*I_2;

    dG_1 = R_Gss_ - hill(X_A,m_G1,h_G1,n_G1).*G_1;
    dG_2 = hill(X_A,m_G1,h_G1,n_G1).*G_1 - hill(X_A,m_G2,h_G2,n_G2).*G_2;
    R_G = hill(X_A,m_G2,h_G2,n_G2).*G_2;

    
    %Mass balances - perfusion
    dI = Q/V_P.*(I_t_in(t) - I) + R_I./V_P;
    dG = Q/V_P.*(G_t_in(t) - G) + R_G./V_P; 
    dg = Q./V_P.*(g_t_in(t) - g); %No glucose generation

    
    dydt = [dI;dG;dX_gB;dX_G;dI_1;dI_2;dX_gA;dX_I;dG_1;dG_2;dg];


end

% Additional functions
function s = R_Gss(X_a,m_a,h_a,n_a)
    %R_Gss represents the steady-state glucagon secretion function
    s = hill(X_a,m_a,h_a,n_a); %mg/min/islet
end 

function s = R_Iss(X_b,m_b,h_b,n_b)
    %R_Iss represents the steady-state insulin secretion function
    s = hill(X_b,m_b,h_b,n_b); %mg/min/islet
end

function s = Y_A(X_gA,X_I,h_IA,n_IA,X_A0,m_g)
    %Y_A represents the net alpha cell signal function
    s = X_gA - hill(X_I,m_g*X_gA+X_A0,h_IA,n_IA) + X_A0;
end
    
function s = Y_B(X_gB,X_G,m_GB,h_GB,n_GB,h_gB,n_gB,X_B0)
    %Y_B represents the net beta cell signal function
    s = X_gB + hill(X_G,m_GB,h_GB,n_GB)*hill(X_gB,1,h_gB,n_gB) + X_B0;
end


function hi = hill(x,m,h,n)

    %Hill Function

    hi = (x >= 0) .* m./((h./x).^n + 1) + (x < 0) .* 0;
    %If the x value is less than 0, the Hill function should still be 0

end