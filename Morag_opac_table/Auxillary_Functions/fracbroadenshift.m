function B = fracbroadenshift(A,broadensize,shiftsize)
%fracbroadenshift broadens and shifts to fractional values, using linear
%interpolation.
% Assumes shiftsize is positive
% broadens in the negative direction

% input must be: [d1,d2,Nnu] table or [1,Nnu] array
%Requires installation of mtimesx

broadensize = abs(broadensize);
int_lo = ceil(shiftsize-broadensize);     %integer portions of broadensize in each direction
fra_lo = int_lo - (shiftsize-broadensize);      %fractional portions of broadensize
int_hi = floor(shiftsize);  
fra_hi = (shiftsize) - int_hi;

if ndims(A)==3
    d1 =size(A,1);
    d2 = size(A,2);
    lA = length(A);

    % B = reshape(A,[d1*d2,lA])*spdiags( [fra_lo*ones(lA,1), ones(lA,int_hi-int_lo+1), fra_hi*ones(lA,1)]/(broadensize+1), int_lo-1 : int_hi+1 , lA, lA );
    B =  mtimesx(reshape(A,[d1*d2,lA]),spdiags( [fra_lo*ones(lA,1), ones(lA,int_hi-int_lo+1), fra_hi*ones(lA,1)]/(broadensize+1), int_lo-1 : int_hi+1 , lA, lA ));
    B(:,1:int_hi+2) = B(:,max(1,int_hi+2)).*ones(1,1,max(int_hi+2)); % Need to change if fracshift is positive
    B = reshape(B,[d1 d2 lA]);
elseif ndims(A)==2
    %% 
    if size(A,1)>1
        error('wrong array size');
    end
    lA = length(A);
    B =  mtimesx(A,spdiags( [fra_lo*ones(lA,1), ones(lA,int_hi-int_lo+1), fra_hi*ones(lA,1)]/(broadensize+1), int_lo-1 : int_hi+1 , lA, lA ));
    if size(A,1)==1
        B(1:int_hi+2) = B(max(1,int_hi+2)).*ones(1,max(1,int_hi+2)); % Need to change if fracshift is positive
    else
        B(:,1:int_hi+2) = B(:,max(1,int_hi+2)).*ones(1,1,max(1,int_hi+2)); % Need to change if fracshift is positive
    end
end
