function [integral_v]=int_nu_dim(nug0,fg5,dim)
% taken from trapz
perm = [dim:max(ndims(fg5),dim) 1:dim-1];
y = permute(fg5,perm);
m = size(y,1);
x = nug0(:);
z = diff(x,1,1).' * y(1:m,:);
siz = size(y); siz(1) = 1;
integral_v = reshape(z,[ones(1,0),siz]);
