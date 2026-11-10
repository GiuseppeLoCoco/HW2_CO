function v = tCG (H_fun,b,delta, max_it)

    it = 0;
    v = zeros(length(b),1);
    r = b - H_fun(v);
    p = r;

    if norm(r, 2) <= min(0.1, norm(b,2)) * (norm(b,2))
        return;
    end
    
    while it<max_it
        Hp = H_fun(p);

     

        if (p'*Hp) <= 0
            ts = roots([norm(p)^2, 2*(v'*p), norm(v)^2 - delta^2]);
            v = v + max(ts) * p;
            return;
        end
        
        alpha = (r' * r) / (p'*Hp);
       
        if norm(v + alpha*p, 2) >= delta
            ts = roots([norm(p)^2, 2*(v'*p), norm(v)^2 - delta^2]);
            v = v + max(ts) * p;
            return;
        end
        
        v  = v + alpha * p;
        r_old = r;
        r = r - alpha * Hp;
    
        if(norm(r,2)<= min(0.1, norm(b,2))*norm(b,2))
            return;
        end
    
        beta = (r'*r)/(r_old'*r_old);
        p = r + beta*p;
        it = it + 1;  
        
    end
end


