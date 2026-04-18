function [Y,X,T,N] = VarStr(y,c,p)
[T,N] = size(y);
Y = y(p+1:T,:);


X = [];

for j = 1:p
    X = [X y(p+1-j:T-j,:)];
end

if c == 1
    X = [ones(size(X,1),1) X];
end

end