function [x, history] = GD_linesearch(x0, fun, grad_fun, alpha_bar, tol, max_it, x_true)
    x = x0;
    
    history.grad_norm = [];
    history.recon_err = [];
    
    rho = 0.5;
    c = 1e-4;
    
    for it = 1:max_it
        fx = fun(x);
        gx = grad_fun(x);
        g_norm = norm(gx, 2);
        
        history.grad_norm(end+1) = g_norm;
        if nargin >= 7 && ~isempty(x_true)
            r_err = min(norm(x - x_true), norm(x + x_true));
            history.recon_err(end+1) = r_err;
        end
        
        if g_norm < tol
            return;
        end
        
        grad_norm_sq = g_norm^2;
        alpha = alpha_bar;
        
        while fx - fun(x - alpha * gx) < c * alpha * grad_norm_sq
            alpha = rho * alpha;
            if alpha < 1e-12
                break;
            end
        end
        
        x = x - alpha * gx;
    end
end