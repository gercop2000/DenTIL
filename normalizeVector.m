function arr=normalizeVector(vect,inverted)
%NORMALIZEVECTOR maps a vector to the range [0,1]

if nargin<2
    inverted=false;
end

if ~inverted % major = black
    arr = ((vect - min(vect)) / (max(vect) - min(vect)));
else
    arr = ((max(vect) - vect) / (max(vect) - min(vect)));
end 

end