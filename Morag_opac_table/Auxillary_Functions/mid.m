function [ vmid ] = mid(v)
%mid function
%vmid=mid(v,dim)
if size(v,1) == 1
    vmid=(v(1:end-1)+v(2:end))/2;
else
    vmid=(v(1:end-1,:)+v(2:end,:))/2;
end

