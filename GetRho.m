function Rho = GetRho (E, h, a)
    
    dRho = 0.05;
    
    r_max = 1;
    
    r = 0;
    Rho(1) = r;
    while ((r+dRho) <= r_max)
        r = r+dRho;
        
        if sqrVr(E, h, r, a) >= 0
            Rho = [Rho r];
        end
    end
    
    %Удаляем 0
    Rho(1) = [];
        
    function result = sqrVr(E, h, r, a)
        result = (2*E + 2*F(r, a) - h^2/r^2);

    function result = F(r, a)
        result = (a^2+r^2)^(-1/2);
            
        