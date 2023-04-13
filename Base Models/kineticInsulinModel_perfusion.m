function dydt = kineticInsulinModel_perfusion(t,y,p,g_t_in)

% kineticInsulinModel_perfusion runs the simplified kinetic-only insulin
% secretion model in a "perfusion" setting: with in and out flows.
% It takes in a scalar time t, a vector of the current state of the system 
% y, a vector of parameters p, and a function handle g_t_in that represents
% the glucose in concentration as a function of time (can only take in 
% time).
% It returns the derivatives of beta cell glucose signal, insulin pool 1, 
% insulin pool 2, glucose concentration, and insulin concentration at the 
% state described by the function inputs

    %Store parameters
    p = num2cell(p);
    [gba, ...
     k_gB, ...
     X_B0, ...
     m_I, h_I, n_I, ...
     m_I1, h_I1, n_I1, m_I2, h_I2, n_I2, ...
     Q, V_P] = p{:};

    %Obtain current system values
    X_gB = y(1); %Glucose signal in beta cells
    I_1 = y(2); %Mass of insulin in pool 1
    I_2 = y(3); %Mass of insulin in pool 2

    g = y(4); %Glucose concentration
    I = y(5); %Inuslin concentration



    %Signal transduction - Assumes that X_G ~ 0, so glucose is the only
    %important signal (Insulin secretion is high, so glucagon is heavily
    %supressed)

    dX_gB = k_gB.*(g./gba - X_gB);

    %Net signals
    X_B = Y_B(X_gB,X_B0);

    %Steady-state secretion
    S_Iss_ = S_Iss(X_B,m_I,h_I,n_I);

    %Transient secretion
    dI_1 = S_Iss_ - hill(X_B,m_I1,h_I1,n_I1).*I_1;
    dI_2 = hill(X_B,m_I1,h_I1,n_I1).*I_1 - hill(X_B,m_I2,h_I2,n_I2).*I_2;
    S_I = hill(X_B,m_I2,h_I2,n_I2).*I_2;

    
    %Mass balances - perfusion
    %No insulin is flowing in, glucagon secretion is assumed to be
    %surpressed and none flows in, so dG = 0 and isn't needed
    dI = Q/V_P.*(0 - I) + S_I./V_P;
    dg = Q./V_P.*(g_t_in(t) - g);

    
    dydt = [dX_gB;dI_1;dI_2;dg;dI];


end

% Additional functions
function s = S_Iss(X_b,m_b,h_b,n_b)
    %S_Iss represents the steady-state insulin secretion function
    s = hill(X_b,m_b,h_b,n_b); %mg/min/islet
end
    
function s = Y_B(X_gB,X_B0)
    %Y_B represents the net beta cell signal function
    s = X_gB + X_B0; %Simplification due to low glucagon levels
end


function hi = hill(x,m,h,n)

    %Hill Function

    hi = (x >= 0) .* m./((h./x).^n + 1) + (x < 0) .* 0;
    %If the x value is less than 0, the Hill function should still be 0

end