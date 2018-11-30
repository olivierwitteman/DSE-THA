function [D_p, w_ee, l_ee, h_ee] = engine_dim_func(P, N)

w_eng = 0.17*(P/(1000*N)).^(0.3);
h_eng = 0.30*(P/(1000*N)).^(0.1);
l_eng = 0.06*(P/(1000*N)).^(0.55);

D_p = 0.55*(P/(1000*N)).^(1/4);

w_ee = 1.2*w_eng;
l_ee = l_eng + 0.1*w_eng;
h_ee = h_eng + 0.2*w_eng;
disp("Engine dimensions are: " + string(w_ee) + " x " + string(l_ee)+ " x " + string(h_ee))
disp("Propeller diameter: " + string(D_p))
end