function [x, y, Vx, Vy] = getIC_full(r, w, E, h, a)
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here
    
    Vw = h/r;
    Vr = sqrt(2*E - Vw^2 + 2*F(r, a));
    
    [x, y] = pol2cart (w, r);
    
    %???
    %Vx = Vr*cos(w) - Vw*sin(w);
    %Vy = Vr*sin(w) + Vw*cos(w);
    
    %modified on 07/01/2009
    Vx = Vr*cos(w) + Vw*sin(w);
    Vy = Vr*sin(w) + Vw*cos(w);
    
    function result = F(r, a)
        result = (a^2+r^2)^(-1/2);