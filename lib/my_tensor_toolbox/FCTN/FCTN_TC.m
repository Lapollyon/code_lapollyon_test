function [X, G] = FCTN_TC(F,Omega,opts)
if isfield(opts, 'rho');         rho   = opts.rho;              end
if isfield(opts, 'R');           R     = opts.R;                end
if isfield(opts, 'maxit');       maxit = opts.maxit;            end
Omega = Omega>0;

X = F;
N = ndims(X); 
Nway = size(X);

tempdim = diag(Nway)+R+R';
G = cell(1,N);
for i = 1:N
    G{i} = rand(tempdim(i,:));
end

for k = 1:maxit
    % Update G 
    Xold = X;
    for i = 1:N
        Xi = my_Unfold(X,Nway,i);
        Gi = my_Unfold(G{i},tempdim(i,:),i);
        Girest = tnreshape(tnprod_rest(G,i),N,i);
        tempC = Xi*Girest'+rho*Gi;
        tempA = Girest*Girest'+rho*eye(size(Gi,2));
        G{i}  = my_Fold(tempC*pinv(tempA),tempdim(i,:),i);
    end
    % Update X 
    X = (tnprod(G)+rho*Xold)/(1+rho);
    X(Omega) = F(Omega);
end
