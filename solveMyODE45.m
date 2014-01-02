function [t_out, Y_out] = solveMyODE45 (fhandle, t, Y, gamma, a)

    [t_out, Y_out] = ode45(fhandle, t, Y, [], [gamma, a]);

end