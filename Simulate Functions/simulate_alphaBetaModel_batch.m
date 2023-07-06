function [t,y] = simulate_alphaBetaModel_batch(params,g_t,I_0,G_0,t)

%simulate_alphaBetaModel_batch runs a simulation of complete alpha-cell and
%beta-cell model in a batch setting. It takes in the parameters for the 
%model as a vector, the glucose concentration trajectory of the system as a
%function, the initial insulin concentration as a scalar, the initial 
%glucagon concentration as a scalar, and the time values that results are 
%desired for as a vector. 
%It returns the time values that have results as a vector and the results 
%of the simulation, including net signals and secretion rates that were 
%calculated after the simulation, as an array.

    %Unpack parameters
    params = num2cell(params);

    [gba_, Gba_, Iba_, ...
     ~, ~, ~, ~, ...
     m_GB_, h_GB_, n_GB_, h_gB_, n_gB_, X_B0_, ...
     h_IA_, n_IA_, X_A0_, m_g_, ...
     m_I_, h_I_, n_I_, ...
     m_G_, h_G_, n_G_, ...
     m_I1_, h_I1_, n_I1_, m_I2_, h_I2_, n_I2_, ...
     m_G1_, h_G1_, n_G1_, m_G2_, h_G2_, n_G2_, ...
     ~] = params{:};

    params = cell2mat(params);

    %Initial glucose concentrations (insulin and glucagon are inputs)
    g_0 = g_t(0);

    %Initial signal intensities
    X_gB_0 = g_0./gba_;
    X_G_0 = G_0./Gba_;
    X_gA_0 = g_0./gba_;
    X_I_0 = I_0./Iba_;

    %Initial net signals
    X_B_0 = Y_B(X_gB_0,X_G_0,m_GB_,h_GB_,n_GB_,h_gB_,n_gB_,X_B0_);
    X_A_0 = Y_A(X_gA_0,X_I_0,h_IA_,n_IA_,X_A0_,m_g_);

    %Initial steady-state secretion rates
    R_Iss_0 = R_Iss(X_B_0,m_I_,h_I_,n_I_);
    R_Gss_0 = R_Gss(X_A_0,m_G_,h_G_,n_G_);

    %Initial pool masses
    I_1_0 = R_Iss_0./hill(X_B_0,m_I1_,h_I1_,n_I1_);
    I_2_0 = hill(X_B_0,m_I1_,h_I1_,n_I1_).*I_1_0./hill(X_B_0,m_I2_,h_I2_,n_I2_);

    G_1_0 = R_Gss_0./hill(X_A_0,m_G1_,h_G1_,n_G1_);
    G_2_0 = hill(X_A_0,m_G1_,h_G1_,n_G1_).*G_1_0./hill(X_A_0,m_G2_,h_G2_,n_G2_);

    %Store initial condititions and time values for integration
    U0_ = [I_0;G_0; ...
           X_gB_0;X_G_0;I_1_0;I_2_0; ...
           X_gA_0;X_I_0;G_1_0;G_2_0];
    t0 = min(t);
    tMax = max(t);

    %Solve system of odes and evaluate at t
    sol = ode23s(@(t,y) alphaBetaModel_batch(t,y,params,g_t) ,[t0 tMax], U0_);
    y = deval(sol,t)';

    %Store signals and second pools to calculate net signals and secretion
    %rates
    X_gB = y(:,3);
    X_G = y(:,4);
    I_2 = y(:,6);

    X_gA = y(:,7);
    X_I = y(:,8);
    G_2 = y(:,10);


    %Calculate X_B, X_A, R_I, and R_G

    %Make placeholder arrays
    X_B = zeros(length(t),1);
    R_I = zeros(length(t),1);
    X_A = zeros(length(t),1);
    R_G = zeros(length(t),1);

    for i = 1:length(t) %Cycle over time values to calculate signals and 
                        %secretion at each time

        X_B(i) = Y_B(X_gB(i),X_G(i),m_GB_,h_GB_,n_GB_,h_gB_,n_gB_,X_B0_);
        R_I(i) = hill(X_B(i),m_I2_,h_I2_,n_I2_).*I_2(i);

        X_A(i) = Y_A(X_gA(i),X_I(i),h_IA_,n_IA_,X_A0_,m_g_);
        R_G(i) = hill(X_A(i),m_G2_,h_G2_,n_G2_).*G_2(i);
    end

    %Store the calculated net signals and secretion with the results of the
    %system of odes to return from this function
    y = [y X_B R_I X_A R_G];


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