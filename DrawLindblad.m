function DrawLindblad( axes_Handle, a )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    hold all;
    
    f1 = sprintf('sqrt((r^4)*((%0.5g)^2+r^2)^(-3/2))', a); %I^(1/2)
    f11 = sprintf('-sqrt((r^4)*((%0.5g)^2+r^2)^(-3/2))', a);
    f2 = sprintf('(r^2/2)*((%0.5g)^2+r^2)^(-3/2)-((%0.5g)^2+r^2)^(-1/2)', a, a); %H
    
    ezplot(axes_Handle, f1, f2);
    ezplot(axes_Handle, f11, f2);
    
    hold off;

end

