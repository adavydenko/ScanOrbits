function dy = my3dODE_full(t, y, varargin)
%   y = {x, y, z, dx, dy, dz}

    gamma = varargin{1}(1);
    a = varargin{1}(2);
    k_r = varargin{1}(3);
    k_z = varargin{1}(4);

    dy = zeros(6, 1);
    
    
    Fx = -(a^2+y(1)^2+y(2)^2+y(3)^2)^(-3/2)*y(1);
    Fy = -(a^2+y(1)^2+y(2)^2+y(3)^2)^(-3/2)*y(2);
    Fz = -(a^2+y(1)^2+y(2)^2+y(3)^2)^(-3/2)*y(3);
    
    dy(1) = y(4);
    dy(2) = y(5);
    dy(3) = y(6);
    dy(4) = +gamma * y(5) + Fx + y(1);
    dy(5) = -gamma * y(4) + Fy;
    dy(6) = Fz - (k_z/k_r)^2 * y(3);