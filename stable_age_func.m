function F = stable_age_func(p)
    global P
    global a
    da = a(2) - a(1);
    F = da.*trapz(P.gH.*exp(-p*a-P.muH_int)) - 1;
end