function  [Y]  =  Im2Patch4D(X, par)
% get full band patches
patsize     = par.patsize;
if isfield(par,'Pstep')
    step   = par.Pstep;
else
    step   = 1;
end
TotalPatNum =  (ceil((size(X,1)-patsize)/step)+1)*(ceil((size(X,2)-patsize)/step)+1);                    %Total Patch Number in the image
Y           =   zeros(patsize, patsize, size(X,3), TotalPatNum);                                             %Patches in the original noisy image

for i  = 1:patsize
    for j  = 1:patsize
        x = [i:step:size(X,1)-patsize+i,size(X,1)-patsize+i]; x= unique(x);
        y = [j:step:size(X,2)-patsize+j,size(X,2)-patsize+j]; y= unique(y);
        tempPatch     =  X(x,y,:);
        Y(i,j,:,:)    =  Unfold(tempPatch, size(tempPatch), 3);
    end
end         %Estimated Local Noise Level