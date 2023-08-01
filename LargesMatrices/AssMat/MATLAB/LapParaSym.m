function A = LapParaSym(Mat11,Mat12,Mat21,Mat22)
% -Parallel_Laplacian Operator

global Dx Dy Nxy
Nx=Nxy(1);
Ny=Nxy(2);

% 0.0625 = 1/16
SqDx = 0.0625/(Dx*Dx);
SqDy = 0.0625/(Dy*Dy);
DxDy = 0.0625/(Dy*Dx);  

V=zeros(3,(Nx-2)*Ny);
iCpt = 1;
Buffer.ij  = zeros(3,3);
Buffer.Val = zeros(3,3);
Val = zeros(2,1);

iCpt = 1;
% Interior Nodes (for I all Nodes for J)
for i=2:(Nx-1)    
    for j=1:(Ny)
        switch j
            case 1
                %% -(qy)_{i,j+1/2} =  -(qy)_{i+1/2,j+1/2} - (qy)_{i-1/2,j+1/2}
                for iMean = 0:1 %0
                    iM = i - iMean;
                    Val(1) = -( Mat21(iM+1,j+1)+Mat21(iM,j+1) ...
                        + Mat21(iM+1,j)  +Mat21(iM,j) ) * DxDy;
                    Val(2) = -( Mat22(iM+1,j+1)+Mat22(iM,j+1) ...
                        + Mat22(iM+1,j)  +Mat22(iM,j) ) * SqDy;
                    Buffer = Qy(Buffer,Val,[1,0]-iMean,[1,0]);
                end
            case Ny
                %% (qy)_{i,j-1/2} = (qy)_{i+1/2,j-1/2} + (qy)_{i-1/2,j-1/2}
                for iMean = 0:1 %1
                    iM = i - iMean;
                    Val(1) = ( Mat21(iM+1,j-1)+Mat21(iM,j-1) ...
                        + Mat21(iM+1,j)  +Mat21(iM,j) ) * DxDy;
                    Val(2) = ( Mat22(iM+1,j-1)+Mat22(iM,j-1) ...
                        + Mat22(iM+1,j)  +Mat22(iM,j) ) * SqDy;
                    Buffer = Qy(Buffer,Val,[1,0]-iMean,[0,-1]);
                end
            otherwise
                % Flux in the X Direction
                %%  (qx)_{i+1/2,j} =  (qx)_{i+1/2,j+1/2} + (qx)_{i+1/2,j-1/2}
                for jMean = 0:1
                    jM = j-jMean;
                    Val(1) = ( Mat11(i+1,jM+1)+Mat11(i,jM+1) ...
                        + Mat11(i+1,jM) + Mat11(i,jM) ) * SqDx;
                    Val(2) = ( Mat12(i+1,jM+1)+Mat12(i,jM+1) ...
                        + Mat12(i+1,jM)  +Mat12(i,jM) ) * DxDy;
                    Buffer = Qx(Buffer,Val,[1,0],[1,0]-jMean);
                    %% -(qx)_{i-1/2,j} =  -(qx)_{i-1/2,j+1/2}-(qx)_{i-1/2,j-1/2}
                    Val(1) = -( Mat11(i-1,jM+1)+Mat11(i,jM+1) ...
                        + Mat11(i-1,jM)  +Mat11(i,jM) ) * SqDx;
                    Val(2) = -( Mat12(i-1,jM+1)+Mat12(i,jM+1) ...
                        + Mat12(i-1,jM)  +Mat12(i,jM) ) * DxDy;
                    Buffer = Qx(Buffer,Val,[0,-1],[1,0]-jMean);
                end
                %%  (qy)_{i,j+1/2} =  (qy)_{i+1/2,j+1/2} + (qy)_{i-1/2,j+1/2}
                for iMean = 0:1
                    iM = i - iMean;
                    Val(1) = ( Mat21(iM+1,j+1)+Mat21(iM,j+1) ...
                        + Mat21(iM+1,j)  +Mat21(iM,j) ) * DxDy;
                    Val(2) = ( Mat22(iM+1,j+1)+Mat22(iM,j+1) ...
                        + Mat22(iM+1,j)  +Mat22(iM,j) ) * SqDy;
                    Buffer = Qy(Buffer,Val,[1,0]-iMean,[1,0]);
                    %% -(qy)_{i,j-1/2} = -(qy)_{i+1/2,j-1/2} - (qy)_{i-1/2,j-1/2}
                    Val(1) = -( Mat21(iM+1,j-1)+Mat21(iM,j-1) ...
                        + Mat21(iM+1,j)  +Mat21(iM,j) ) * DxDy;
                    Val(2) = - ( Mat22(iM+1,j-1)+Mat22(iM,j-1) ...
                        + Mat22(iM+1,j)  +Mat22(iM,j) ) * SqDy;
                    Buffer = Qy(Buffer,Val,[1,0]-iMean,[0,-1]);
                end
        end
        Numl = lexico(i,j);
        for j0 = -1:1
            jb = j0+2;
            for i0 = -1:1
                ib = i0+2;
                if Buffer.ij(ib,jb) == 1
                    % Store current line
                    Numc = lexico(i+i0,j+j0);
                    V(1,iCpt) = Numl;
                    V(2,iCpt) = Numc;
                    V(3,iCpt) = Buffer.Val(ib,jb);
                    iCpt = iCpt + 1 ;
                    %Reset Buffer
                    Buffer.ij(ib,jb) = 0 ;
                    Buffer.Val(ib,jb) = 0 ;
                end
            end
        end
        
    end
end
A=sparse(V(1,:),V(2,:),V(3,:),Nx*Ny,Nx*Ny);
    


function Buffer=Qx(Buffer,Val,I,J)
    % {i+1,j+1/2}
    Buffer = StackBuffer(Buffer,I(1),J(1),-Val(1));
    Buffer = StackBuffer(Buffer,I(1),J(2),-Val(1));
    % {i,j+1/2}
    Buffer = StackBuffer(Buffer,I(2),J(1), Val(1));
    Buffer = StackBuffer(Buffer,I(2),J(2), Val(1));        
    % {i+1/2,j+1}
    Buffer = StackBuffer(Buffer,I(1),J(1),-Val(2));
    Buffer = StackBuffer(Buffer,I(2),J(1),-Val(2));
    % {i+1/2,j  }
    Buffer = StackBuffer(Buffer,I(1),J(2), Val(2));
    Buffer = StackBuffer(Buffer,I(2),J(2), Val(2)); 
end

function Buffer=Qy(Buffer,Val,I,J)
    % {i+1,j+1/2}
    Buffer = StackBuffer(Buffer,I(1),J(1),-Val(1));
    Buffer = StackBuffer(Buffer,I(1),J(2),-Val(1));
    % {i,j+1/2}
    Buffer = StackBuffer(Buffer,I(2),J(1), Val(1));
    Buffer = StackBuffer(Buffer,I(2),J(2), Val(1));
            
    % {i+1/2,j+1}
    Buffer = StackBuffer(Buffer,I(1),J(1),-Val(2));
    Buffer = StackBuffer(Buffer,I(2),J(1),-Val(2));
    % {i+1/2,j  }
    Buffer = StackBuffer(Buffer,I(1),J(2), Val(2));
    Buffer = StackBuffer(Buffer,I(2),J(2), Val(2));    
end
end

