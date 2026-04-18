function [C,ex_out] = companion(A,p,N,ex)
    arguments
        A
        p
        N
        ex = 0;
    end
    
 
    AR_coeffs = A(2:1+(N*p), :); 
    M = AR_coeffs.';             
    

    if p > 1

        C = [M; 
             eye(N*(p-1)), zeros(N*(p-1), N)];
    else
        C = M; 
    end
    
    ex_out = ex; 
end



