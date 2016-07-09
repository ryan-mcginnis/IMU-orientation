function q = dcm2quatern(R)
%Function to extract a quaternion from a DCM of size (3,3,n).  
%Method from Bar-Itzhack 2000
    numR = size(R,3);
    q = zeros(numR, 4);
    K = zeros(4,4);
    for i = 1:numR
        K(1,1) = (1/3) * (R(1,1,i) - R(2,2,i) - R(3,3,i));
        K(1,2) = (1/3) * (R(2,1,i) + R(1,2,i));
        K(1,3) = (1/3) * (R(3,1,i) + R(1,3,i));
        K(1,4) = (1/3) * (R(2,3,i) - R(3,2,i));
        K(2,1) = (1/3) * (R(2,1,i) + R(1,2,i));
        K(2,2) = (1/3) * (R(2,2,i) - R(1,1,i) - R(3,3,i));
        K(2,3) = (1/3) * (R(3,2,i) + R(2,3,i));
        K(2,4) = (1/3) * (R(3,1,i) - R(1,3,i));
        K(3,1) = (1/3) * (R(3,1,i) + R(1,3,i));
        K(3,2) = (1/3) * (R(3,2,i) + R(2,3,i));
        K(3,3) = (1/3) * (R(3,3,i) - R(1,1,i) - R(2,2,i));
        K(3,4) = (1/3) * (R(1,2,i) - R(2,1,i));
        K(4,1) = (1/3) * (R(2,3,i) - R(3,2,i));
        K(4,2) = (1/3) * (R(3,1,i) - R(1,3,i));
        K(4,3) = (1/3) * (R(1,2,i) - R(2,1,i));
        K(4,4) = (1/3) * (R(1,1,i) + R(2,2,i) + R(3,3,i));
        [V,~] = eig(K);
        q(i,:) = V(:,4)';
        q(i,:) = [q(i,4) q(i,1) q(i,2) q(i,3)] ./ norm(q(i,:));
        
        %Check for continuity
        if numR>1 && i>1
            dq = sum(abs(q(i,:) - q(i-1,:)));
            dqp = sum(abs(-q(i,:) - q(i-1,:)));
            
            if dq > dqp
                q(i,:) = -q(i,:);
            end
        end
    end
end