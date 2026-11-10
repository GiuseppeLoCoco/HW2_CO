function x = TR_tCG(x0, delta_0, fun, H_fun, b_fun, rho_prime, delta_max, max_it)
    it = 0;
    delta = delta_0;
    x = x0;
    
    while (it < max_it)
        b = b_fun(x); 
        H_at_x = @(u) H_fun(x, u);
        
        uk = tCG(x0, H_at_x, b, delta, max_it);
        x_star = x + uk;
        
        Huk = H_at_x(uk); 
        delta_m = -(0.5 * (uk' * Huk) + b' * uk); 
        
        rho_k = (fun(x) - fun(x_star)) / delta_m;
        
        if (rho_k > rho_prime)
            x = x_star;
        end
        
        if (rho_k < 1/4)
            delta = delta / 4;
        elseif (rho_k > 3/4 && norm(uk) == delta)
            delta = min(2 * delta, delta_max);
        end
        
        it = it + 1;
    end
end