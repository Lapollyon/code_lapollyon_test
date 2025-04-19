function  [E_Img, Weight]  =  Patch4D2Im(ImPat, WPat, par, sizeV)
patsize      = par.patsize;
if isfield(par,'Pstep')
    step   = par.Pstep;
else
    step   = 1;
end
TempR        =   ceil((sizeV(1)-patsize)/step)+1;
TempC        =   ceil((sizeV(2)-patsize)/step)+1;

E_Img  	=  zeros(sizeV);
W_Img 	=  zeros(sizeV);

for i  = 1:patsize
    for j  = 1:patsize
        x = [i:step:sizeV(1)-patsize+i,sizeV(1)-patsize+i]; x= unique(x);
        y = [j:step:sizeV(2)-patsize+j,sizeV(2)-patsize+j]; y= unique(y);
        E_Img(x,y,:)  =  E_Img(x,y,:) + Fold(ImPat(i,j,:,:),[TempR TempC sizeV(3)], 3);
        W_Img(x,y,:)  =  W_Img(x,y,:) + Fold(repmat(WPat(i,j,:),sizeV(3),1), [TempR TempC sizeV(3)], 3);
    end
end
E_Img  =  E_Img./(W_Img+eps);
Weight =  1./(W_Img);

