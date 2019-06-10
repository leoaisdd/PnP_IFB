function y = afun(x,transp_flag,h,dim)

rows = dim(1);
cols = dim(2);
if strcmp(transp_flag,'transp')         % y = A'*x
    x = reshape(x,[rows,cols]);
    y = imfilter(x,rot90(h,2),'circular');
    y = y(:);
elseif strcmp(transp_flag,'notransp')   % y = A*x
    x = reshape(x,[rows,cols]);
    y = imfilter(x,h,'circular');
    y = y(:);
end
end