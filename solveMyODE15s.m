function [t_out, Y_out] = solveMyODE15s (t, Y, gamma, r0)

    [t_out, Y_out] = ode15s(@myODE, t, Y, [], [gamma, r0]);

end