function drawHill (r0, C)
%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here
f='';
f = sprintf('2*(1+(x^2+y^2)/(%d^2))^(-1/2)+x^2-%d', r0, C);

ezplot(f);

end