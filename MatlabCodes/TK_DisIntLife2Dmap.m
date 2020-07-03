function Out = TK_DisIntLife2Dmap(Int, Life, Z2Dmap, SmoothFactor, RemoveLife0orNot)
% All points with lifetime of 0 become Z = 0
% RemoveLife0orNot = 1  =>  Remove all points with lifetime of 0

A = Z2Dmap ;
if (RemoveLife0orNot == 1)
    ImaxLife = length(Life) ;
    B = Life(2:ImaxLife) - Life(1:ImaxLife-1) ;
    LimitLife = mean(B) .* 10^-2 ;
    I = 0 ;
    while I < ImaxLife
        I = I + 1 ;
        if (Life(I) >= -1.*LimitLife) && (Life(I) <= LimitLife)
            A(I,:) = 0 ;
        end
    end
end

figure
%imgaussfilt(IntLife2Dcor,n) => 2D smoothing by factor of n
image(Int, Life, imgaussfilt(A,SmoothFactor), 'CDataMapping','scaled') ;
axis('xy') ;        % if axis propety needs, enter " [S1,S2,S3] = axis('state') "
load('mycmap.dat');
colormap(mycmap);

global g2Dmax ;
g2Dmax = max(max(imgaussfilt(A,SmoothFactor))) ;

clear A B LimitLife;
end