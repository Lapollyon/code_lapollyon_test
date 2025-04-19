function [Emsi] = NL_FCTN_TC(Nmsi, Omsi, Omega, par, parP)
% MSI denoising based on FCTN
%==========================================================================
if nargin < 4
    OutPSNR = 0;
else
    OutPSNR = 1;
end


sizeData        = size(Nmsi);
Npatch          = Im2Patch4D(Nmsi, par);
sizePatch       = size(Npatch);
Emsi            = Nmsi;
[Sel_arr]       = nonLocal_arr(sizeData, par); % PreCompute the all the patch index in the searching window
L               = length(Sel_arr);
%% main loop
tic
    for iter = 1 : par.deNoisingIter
        Curmsi      	= Emsi;
        Curpatch        = Im2Patch4D(Curmsi, par); % image to FBPs...
        Omegapatch      = Im2Patch4D(Omega, par);
        
        % %    block matching to find samilar FBP goups
        unfoldPatch     = Unfold(Curpatch,sizePatch,4)';
        patchXpatch     = sum(unfoldPatch.^2,1);
        index     =  zeros(par.patnum,L); 
        sizePart  =  250; % This can be change according to memory
        numPart   =  floor(L/sizePart)+1;
        fprintf('Block matching of iter %f has been done          ', iter);
        fprintf('\n')
            for i = 1:numPart
                tempInd = (i-1)*sizePart+1:min(L,i*sizePart);
                distenMat         = bsxfun(@plus, patchXpatch(Sel_arr(tempInd)), patchXpatch')-2*(unfoldPatch')*unfoldPatch(:,Sel_arr(tempInd));
                [~,tempindex]     = sort(distenMat);
                index(:,tempInd)  = tempindex(1:par.patnum,:);
                if i/numPart<0.1
                    fprintf('\b\b\b\b\b\b\b %2.2f%% ', i/numPart*100);
                else
                    fprintf('\b\b\b\b\b\b\b\b %2.2f%% ', i/numPart*100);
                end
            end
            fprintf('\n')
            clear patchXpatch distenMat unfoldPatch Emsi Curmsi tempInd ;
            Epatch           = zeros(sizePatch);
            W                = zeros(sizePatch(1),sizePatch(2),sizePatch(4),'single');
            fprintf('NL-FCTN of iter %f has been done          ', iter);
            fprintf('\n')
            sizePart  =  100;
            numPart   =  floor(L/sizePart)+1;
        for i = 1:numPart
            PattInd = (i-1)*sizePart+1:min(L,i*sizePart);
            tempInd = index(:,PattInd);
            sizeInd = size(tempInd);
            tempPatch = Curpatch(:,:,:,tempInd(:));
            tempOmegapatch = Omegapatch(:,:,:,tempInd(:));
            tempPatch = reshape(tempPatch, [sizePatch(1:3), sizeInd]);
            tempOmegapatch = reshape(tempOmegapatch, [sizePatch(1:3), sizeInd]);
            parfor j = 1:sizeInd(2)
                [tempPatch(:,:,:,:,j),~] = FCTN_TC(tempPatch(:,:,:,:,j), tempOmegapatch(:,:,:,:,j), parP); % Perform FCTN-based tensor recovery on each FBP goup
            end
            for j = 1:sizeInd(2)
                Epatch(:,:,:,tempInd(:,j))  = Epatch(:,:,:,tempInd(:,j)) + tempPatch(:,:,:,:,j);
                W(:,:,tempInd(:,j))         = W(:,:,tempInd(:,j))+ones(size(tempPatch(:,:,:,:,j),1),size(tempPatch(:,:,:,:,j),2),size(tempPatch(:,:,:,:,j),4));
            end
            if i/numPart<0.1
                fprintf('\b\b\b\b\b\b\b %2.2f%% ', i/numPart*100);
            else
                fprintf('\b\b\b\b\b\b\b\b %2.2f%% ', i/numPart*100);
            end
        end
        clear Curpatch;
        fprintf('\n')
        
    
    clear tempPatch
    time = toc;
    [Emsi, ~]  =  Patch4D2Im(Epatch,W,par,sizeData); % recconstruct the estimated MSI by aggregating all reconstructed FBP goups.
    clear Epatch;
    
    if OutPSNR
        psnr       =  my_PSNR(Emsi*255,Omsi*255);
        disp(['Iter: ' num2str(iter),' , current PSNR = ' num2str(psnr), ',  already cost time: ', num2str(time)]);
    else
        disp(['Iter: ' num2str(iter),'   done,  already cost time: ', num2str(time)]);
    end

    end
end


function  [SelfIndex_arr]  =  nonLocal_arr(sizeD, par)
% -SelfIndex_arr is the index of keypatches in the total patch index array
TempR         =   sizeD(1)-par.patsize+1;
TempC         =   sizeD(2)-par.patsize+1;
R_GridIdx	  =   [1:par.step:TempR,TempR];
R_GridIdx     =   unique(R_GridIdx);
C_GridIdx	  =   [1:par.step:TempC,TempC];
C_GridIdx     =   unique(C_GridIdx);

temp          = 1:TempR*TempC;
temp          = reshape(temp,TempR,TempC);
SelfIndex_arr = temp(R_GridIdx,C_GridIdx);
SelfIndex_arr = SelfIndex_arr(:)';
end