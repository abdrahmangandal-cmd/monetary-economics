%OLS
function [beta,res,T,N,X] = OLS(data,p,ex,d)
    arguments
        data
        p
        ex = 0;
        d = 0; 
    end


[Y,X,T,N] = VarStr(data,1,p);


if isscalar(ex) && ex==0
    X = X;
    if isscalar(d) && d == 0
        X = X;
    else
        X = [X X(:,2:end).*d(p+1:end,:)];
    end
else
    if isscalar(d) && d == 0
        X = [X ex(p+1:end,:)];
    else
        X = [X X(:,2:end).*d(p+1:end,:) d(p+1:end,:) ex(p+1:end,:) ex(p+1:end,:).*d(p+1:end,:)];

    end    
end

beta = (X.'*X) \ (X.'*Y);
res = Y-X*beta;
 



end


%gives back the constant (1+Np)xN