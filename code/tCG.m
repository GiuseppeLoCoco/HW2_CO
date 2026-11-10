function v = tCG (v0,H_fun,b,delta, max_it)

    it = 0;
    v = v0;
    r = b - H_fun(v);
    p = r;
    
    while it<max_it
        Hp = H_fun(p);
        alpha = r' * r/(p'*Hp);
    
        if((p'*Hp)<=0 || norm((v + alpha*p),2)>=delta)
            ts = roots([norm(p)^2, 2*(v'*p), norm(v)^2-delta^2]);
            t_plus = max(ts);
            v  = v + t_plus * p;   
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


