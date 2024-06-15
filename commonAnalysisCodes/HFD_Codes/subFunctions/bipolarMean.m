function M = bipolarMean(A,changeDimFlag,medianFlag)
% this function is used to get the mean for bipolar electrodes such that the
% final matrix is 3x3x:
% if changeDimFlag, then final matrix is 9x:
 
% it does mean across first dimension
if ~exist('changeDimFlag','var')    changeDimFlag=1;    end
if ~exist('medianFlag','var')    medianFlag=0;    end
z = size(A);

if mod(z(1),3)~=0
    disp('mean operation cannot be performed. Wrong size of the matrix')
    M=A;
    return;
else
    M = zeros([3 z(2:length(z))]);
    z1  = 1: z(1)/3;

    for i = 1:3
        if medianFlag
             M(i,:) = nanmean(A([(z1-1)*3+i],:),1);
        else
            M(i,:) = nanmean(A([(z1-1)*3+i],:),1);
        end
    end

    if changeDimFlag
      M = squeeze(reshape(M,1,9,[]));
    end
end
end
        