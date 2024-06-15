function yinterp = interpData(x,y,xinterp)
yinterp = [];
%y is 3 D dimensional data and is interpolated across 3rd dimension
    if length(size(y))~=3
        disp('wrong dimension of y');
        return;
    end

    for i=1:size(y,1)
        for j=1:size(y,2)
            if ~isnan(x(i,j,:))
                 yinterp(i,j,:) = interp1(squeeze(x(i,j,:)),squeeze(y(i,j,:)),xinterp,'pchip','extrap');
            else
                yinterp(i,j,:) = nan(1,length(xinterp));
            end
        end
    end

end