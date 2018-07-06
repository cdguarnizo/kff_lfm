function D = computeD(alpha, omega, lambda)
    
D = - (alpha + 1i*lambda)./(alpha^2 + omega^2 - lambda.^2 + 1i*2*alpha*lambda);