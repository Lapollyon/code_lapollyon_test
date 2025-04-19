function Out = tnprod_rest(X,rest)
N = length(X);
m1 = zeros(1,N-2); m2 = zeros(1,N-2);
n = 1:N-1; 
if rest <N
   n(rest)=[];
end
for i=1:N-2
    m1(i)=2+(i-1)*N;
    m2(i)=3+(i-1)*N;
end
j=1;
if rest>1
    Out = X{1};
    for i=1:N-1
        if i+1 < rest
            Out = ttprod_mn(Out,X{i+1},m1(1:j),n(1:j));
            m1(2:end)=m1(2:end)-[1:N-3];
            m2(2:end)=m2(2:end)-[1:N-3];
            j=j+1;
        end
        if i+1 > rest
            Out = ttprod_mn(Out,X{i+1},m2(1:j),n(1:j));
            m2(2:end)=m2(2:end)-[1:N-3];
            j=j+1;
        end
    end
end
if rest == 1
    Out = X{2};
    for i=2:N-1
        Out = ttprod_mn(Out,X{i+1},m2(1:i-1),n(1:i-1));
        m2(2:end)=m2(2:end)-[1:N-3];
    end
end