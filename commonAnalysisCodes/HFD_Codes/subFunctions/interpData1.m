%Getting 1D interpolation for a matrix
%incomplete code
function yinterp = interpData1(x,y,xinterp,dimForInterp)
yinterp = [];

dimy = length(size(y));
index = cell(1,dimy);
index{dimForInterp} = ':';
analysisMatrix = setdiff([1:dimy],dimForInterp);
index(analysisMatrix) = {1};

for i = flip(analysisMatrix)
    for j = 1:size(y,i)
        index{i} = j;
        disp(index);
    end
end
%y is 3 D dimensional data and is interpolated across 3rd dimension
    if length(size(y))~=3
        disp('wrong dimension of y');
        return;
    end

    for i=1:size(y,1)
        for j=1:size(y,2)
            if ~isnan(x(i,j,:))
                 yinterp(i,j,:) = interp1(squeeze(x(i,j,:)),squeeze(y(i,j,:)),xinterp,'linear','extrap');
            else
                yinterp(i,j,:) = nan(1,length(xinterp));
            end
        end
    end

end