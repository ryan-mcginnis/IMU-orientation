function qConj = quaternConj(q)
%Function to calculate conjugate of quaternion 
%Inputs:
%1. q - input quaternion (nx4)

%Outputs:
%1. qConj - quaternion conjugate of q

qConj = [q(:,1) -q(:,2) -q(:,3) -q(:,4)];

end