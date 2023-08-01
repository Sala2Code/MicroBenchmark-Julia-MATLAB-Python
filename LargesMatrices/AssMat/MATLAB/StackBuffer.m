function Buffer = StackBuffer(Buffer,i,j,Val)
        Buffer.ij(i+2,j+2) = 1;
        Buffer.Val(i+2,j+2) = Buffer.Val(i+2,j+2) + Val;