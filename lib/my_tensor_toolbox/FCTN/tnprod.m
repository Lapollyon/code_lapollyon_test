function Out = tnprod(X)
N=length(X);
m = 2; n = 1;
Out = X{1};
for i=1:N-1
    Out = ttprod_mn(Out,X{i+1},m,n);
    n=[n,1+i];
    tempm = 2+i*(N-i);
    if i>1
        m(2:end)=m(2:end)-[1:i-1];
    end
    m   = [m, tempm];
end