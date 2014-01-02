function [t_out, Y_out] = solveMyODE113 (t, Y, gamma, r0)

    [t_out, Y_out] = ode113(@myODE, t, Y, [], [gamma, r0]);

end