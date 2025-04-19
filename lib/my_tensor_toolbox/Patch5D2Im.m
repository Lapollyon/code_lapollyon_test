function  [E_Img, Weight]  =  Patch5D2Im(ImPat, WPat, par, sizeV)
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
        tempE_Img = reshape(shiftdim(squeeze(ImPat(i,j,:,:,:)), 2), TempR, TempC, sizeV(3), sizeV(4));
        tempW_Img = reshape(shiftdim(repmat(WPat(i,j,:), [sizeV(3),sizeV(4)]), 2), TempR, TempC, sizeV(3), sizeV(4));
        E_Img(x,y,:,:) =  E_Img(x,y,:,:) + tempE_Img;
        W_Img(x,y,:,:) =  W_Img(x,y,:,:) + tempW_Img;
    end
end
E_Img  =  E_Img./(W_Img+eps);
Weight =  1./(W_Img);

