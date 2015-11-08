function dydt = Ebola2(t, y, r, q,beta,omega,rho)
dydt = [-r*q*(y(2)/2)*y(1)/(q*y(2) + y(1) + y(3)) - beta*y(1); 
    r*q*(y(2)/2)*y(1)/(q*y(2) + y(1) + y(3)) - omega*y(2)/2 - (rho)*y(2)/2;
    (rho)*y(2)/2 + beta*y(1); omega*y(2)/2];