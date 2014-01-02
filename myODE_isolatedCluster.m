function dy = myODE_isolatedCluster(t, y, varargin)
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here
    
    dy = zeros(4, 1);

    dy(1) = y(3);
    dy(2) = y(4);
    dy(3) = -y(1) * (1 + y(1)^2 + y(2)^2)^(-3/2);
    dy(4) = -y(2) * (1 + y(1)^2 + y(2)^2)^(-3/2);