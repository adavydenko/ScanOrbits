function [t_out, Y_out] = solveMy3dODE45 (fhandle, t, Y, gamma, a, k_r, k_z)

    [t_out, Y_out] = ode45(fhandle, t, Y, [], [gamma, a, k_r, k_z]);

end