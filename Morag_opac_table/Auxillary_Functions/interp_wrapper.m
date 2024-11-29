function Vq=interp_wrapper(V,Xq,Mode)
% Uses ScaleTime if user chooses to do so, otherwise uses intep1
if Mode.kappa.use_ScaleTime
    Vq = ScaleTime(V,Xq);
else
    Vq = interp1(V,Xq);
end
end