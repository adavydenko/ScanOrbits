function [t_out, Y_out] = solveMyODE23t (t, Y, gamma, r0)

    [t_out, Y_out] = ode23t(@myODE, t, Y, [], [gamma, r0]);

end