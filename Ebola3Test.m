function dydt = Ebola3Test(t,y,r,q,Q,omega,rho)
beta = sawtooth(-1/2*pi*t, 1);
dydt = [-r*q*(y(2)/2)*y(1)/(q*y(2) + y(1) + y(3)) - Q*(beta +1)*y(1);
    r*q*(y(2)/2)*y(1)/(q*y(2) + y(1) + y(3)) - omega*y(2)/2 - (rho)*y(2)/2;
    (rho)*y(2)/2 + Q*(beta+1)*y(1);];
