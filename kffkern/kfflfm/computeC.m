function C = computeC(alpha, omega, lambda)
    
C = omega./(alpha^2 + omega^2 - lambda.^2 + 1i*2*alpha*lambda);