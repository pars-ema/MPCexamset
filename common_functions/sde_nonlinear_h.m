function z = sde_nonlinear_h(X,p,R)
    z = stoch_nonlinear_g(X,p,R);
    z = z(1:2);
end