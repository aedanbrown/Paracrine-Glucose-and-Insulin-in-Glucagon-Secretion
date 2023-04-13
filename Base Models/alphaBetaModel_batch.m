function dydt = alphaBetaModel_batch(t,y,p,g_t)

% alphabetaModel_batch runs the alpha- and beta-cell secretion model in a
% "batch" setting: no in or out flows.
% It takes in a scalar time t, a vector of the current state of the system 
% y, a vector of parameters p, and a function handle g_t that represents
% the glucose concentration as a function of time (can only take in time)
% It returns the derivatives of insulin concentration, glucagon
% concentration, beta cell glucose signal, beta cell glucagon signal,
% insulin pool 1, insulin pool 2, alpha cell glucagon signal, alpha cell
% insulin signal, glucagon pool 1, and glucagon pool 2 at the state
% described by the function inputs

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
     V_P] = p{:};

    %Obtain current system values
    g = g_t(t); %Glucose concentration
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
    S_Iss_ = S_Iss(X_B,m_I,h_I,n_I);
    S_Gss_ = S_Gss(X_A,m_G,h_G,n_G);

    %Transient secretion
    dI_1 = S_Iss_ - hill(X_B,m_I1,h_I1,n_I1).*I_1;
    dI_2 = hill(X_B,m_I1,h_I1,n_I1).*I_1 - hill(X_B,m_I2,h_I2,n_I2).*I_2;
    S_I = hill(X_B,m_I2,h_I2,n_I2).*I_2;

    dG_1 = S_Gss_ - hill(X_A,m_G1,h_G1,n_G1).*G_1;
    dG_2 = hill(X_A,m_G1,h_G1,n_G1).*G_1 - hill(X_A,m_G2,h_G2,n_G2).*G_2;
    S_G = hill(X_A,m_G2,h_G2,n_G2).*G_2;

    
    %Mass balances - batch
    dI = S_I/V_P;
    dG = S_G/V_P; 

    
    dydt = [dI;dG;dX_gB;dX_G;dI_1;dI_2;dX_gA;dX_I;dG_1;dG_2];


end

% Additional functions
function s = S_Gss(X_a,m_a,h_a,n_a)
    %S_Gss represents the steady-state glucagon secretion function
    s = hill(X_a,m_a,h_a,n_a); %mg/min/islet
end 

function s = S_Iss(X_b,m_b,h_b,n_b)
    %S_Iss represents the steady-state insulin secretion function
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