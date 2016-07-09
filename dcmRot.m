function xp = dcmRot(R,x)
%Function to apply rotation defined by direction cosine matrix
%Inputs:
%1. R - direction cosine matrix (3x3xn) defining rotation to be applied
%2. x - 3-element vector (nx3) to be transformed by quaternion

%Outputs:
%1. xp - transformed vector (nx3)

%Account for case where x=(1x3), R=(3x3xn)
if size(x,1)==1 && size(R,3)~=1
    x = ones(size(R,3),1) * x;
end

%Apply rotation
xp = zeros(size(x));
for row_ind = 1:length(xp)
    xp(row_ind,:) = x(row_ind,:) * R(:,:,row_ind).';
end
end