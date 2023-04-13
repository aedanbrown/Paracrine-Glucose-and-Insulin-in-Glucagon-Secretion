function results = simulate_kineticInsulinModel_perfusion_multi(p,g_in_t_all,t)

%simulate_kineticInsulinModel_perfusion_multi runs multiple simulations of
%the simplified insulin secretion model used for curve fitting in a
%perfusion setting with different glucose in flow rate trajectories. It 
%takes in theparameters for each simulation as a vecotr, the glucose in 
%flow rate trajectories as a cell array of functions, and the time values
%that results are desired as a vector.
%It returns a 3-d matrix of the results for each condition.

    results = zeros(length(t),7,length(g_in_t_all));
    %Time along rows, independent variables along columns, conditions along depths

    for i = 1:length(g_in_t_all)

        [~,y] = simulate_kineticInsulinModel_perfusion(p,g_in_t_all{i},t);
        results(:,:,i) = y;

    end

end