function [t,y] = simulate_kineticInsulinModel_perfusion(params,g_t_in,t)

%simulat_kineticInsulinModel_perfusion runs a simulation of the simplified
%insulin secretion model used for curve fitting in a perfusion setting. It
%takes in the parameters for the model as a vector, the glucose in flow
%rate trajectory as a function, and the time values that results are
%desired at as a vector.
%It returns the time values that results are desired at as a vector and the
%results of the simulation, including net signals and secretion rates that
%were calculated after the simulation, as an array.

    %Unpack parameters
    params = num2cell(params);
    [gba_, ...
     ~, ...
     X_B0_, ...
     m_I_, h_I_, n_I_, ...
     m_I1_, h_I1_, n_I1_, m_I2_, h_I2_, n_I2_, ...
     Q_, ~] = params{:};

    params = cell2mat(params);

    
    %Initial glucose in concentration
    g_in_0 = g_t_in(0);

    %Initial signal intensities
    X_gB_0 = g_in_0./gba_;

    %Initial net signals
    X_B_0 = Y_B(X_gB_0,X_B0_);

    %Initial steady-state secretion rates
    S_Iss_0 = S_Iss(X_B_0,m_I_,h_I_,n_I_);

    %Initial pool masses
    I_1_0 = S_Iss_0./hill(X_B_0,m_I1_,h_I1_,n_I1_);
    I_2_0 = hill(X_B_0,m_I1_,h_I1_,n_I1_).*I_1_0./hill(X_B_0,m_I2_,h_I2_,n_I2_);

    %Initial insulin concentration
    I_0 = S_Iss_0./Q_;

    %Store initial conditions and time values for integration
    U0_ = [X_gB_0;I_1_0;I_2_0;g_in_0;I_0];
    t0 = min(t);
    tMax = max(t);

    %Solve system of odes and evaluate at t
    sol = ode23(@(t,y) kineticInsulinModel_perfusion(t,y,params,g_t_in), [t0 tMax], U0_);
    y = deval(sol,t)';

    %Store signals and second pools to calculate net signals and secretion
    %rates
    X_gB = y(:,1);
    I = y(:,5);

    %Calculate X_B and S_I
    X_B = zeros(length(t),1);

    for i = 1:length(t)
         X_B(i) = Y_B(X_gB(i),X_B0_);
    end
    S_I = Q_.*I; %This is the measured secretion rate, which may vary 
                 %slightly from the true beta-cell secretion rate, which is
                 %hill(X_B(i),m_I2_,h_I2_,n_I2_).*I_2(i);

    y = [y X_B S_I];


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