function spread_mn = SpreadProp(y_m,y_n,q_m,q_n)
    N_m = (q_m*y_m(2)+y_m(1)+y_m(3));
    N_n = (q_n*y_n(2)+y_n(1)+y_n(3));
    spread_mn=(q_m*y_m(2)/2)/(N_m)*min(1,N_m^2/N_n)*(y_n(1)/N_n);
end