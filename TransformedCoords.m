function [psi, d_psi, eta, d_eta] = TransformedCoords(A)
    
    r0 = getR0;
    t0 = getT0;

    psi = A(1)/r0;
    eta = A(3)/r0;
    d_psi = A(2)*t0/r0;
    d_eta = A(4)*t0/r0;
    
    function result = getT0
        result = 1/5;
    end

    function result = getR0
        result = 0.8;
    end

end