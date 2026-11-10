function [x, history] = TR_tCG(x0, delta_0, fun, H_fun, b_fun, rho_prime, delta_max, max_it, x_true)
    it = 0;
    delta = delta_0;
    x = x0;
    
    history.grad_norm = [];
    history.recon_err = [];
    
    while (it < max_it)
        b = b_fun(x); 
        
        g_norm = norm(b, 2);
        history.grad_norm(end+1) = g_norm;
        
        if nargin >= 9 && ~isempty(x_true)
            r_err = min(norm(x - x_true), norm(x + x_true));
            history.recon_err(end+1) = r_err;
        end
        
        if g_norm < 1e-5
            return;
        end
        
        H_at_x = @(u) H_fun(x, u);
        
        uk = tCG(H_at_x, b, delta, max_it);
        
        if norm(uk) < 1e-12
            return;
        end
        
        x_star = x + uk;
        
        Huk = H_at_x(uk); 
        delta_m = -0.5 * (uk' * Huk) + b' * uk; 
        
        if delta_m <= 0
            delta = delta / 4;
            it = it + 1;
            continue;
        end
        
        rho_k = (fun(x) - fun(x_star)) / delta_m;
        
        if (rho_k > rho_prime)
            x = x_star;
        end
        
        if (rho_k < 1/4)
            delta = delta / 4;
        elseif (rho_k > 3/4 && norm(uk) >= 0.9 * delta)
            delta = min(2 * delta, delta_max);
        end
        
        it = it + 1;
    end
end