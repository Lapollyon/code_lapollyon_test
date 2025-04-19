function Out = ttprod_mn(X,Y,m,n)
Lx = size(X);      Ly = size(Y);   
Nx = ndims(X);     Ny = ndims(Y);
indexx = 1:Nx;     indexy = 1:Ny;
indexx(m) = [];    indexy(n) = [];

tempX = permute(X,[indexx,m]);  tempXX=reshape(tempX,prod(Lx(indexx)),prod(Lx(m)));
tempY = permute(Y,[n,indexy]);  tempYY=reshape(tempY,prod(Ly(n)),prod(Ly(indexy)));
tempOut = tempXX*tempYY;
Out     = reshape(tempOut,[Lx(indexx),Ly(indexy)]);
