function xp = quaternRot(q,x)
%Function to apply quaternion rotation 
%Inputs:
%1. q - quaternion (nx4) defining rotation to be applied
%2. x - 3-element vector (nx3) to be transformed by quaternion

%Outputs:
%1. xp - transformed vector (nx3)

%Pad x with column of zeros (quaternion format)
x = [zeros(size(x,1),1), x];

%Account for case where x=(1x3), q=(nx4)
if size(x,1)==1 && size(q,1)~=1
    x = ones(size(q,1),1) * x;
end

%Apply rotation
xt = quaternProd(q, quaternProd(x, quaternConj(q)));

%Extract rotated vector
xp = xt(:,2:4);
end