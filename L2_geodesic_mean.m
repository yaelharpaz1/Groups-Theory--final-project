function S = L2_geodesic_mean(R)

% Input:   R, a 3x3xn matrix, represents a set of rotation matrices in SO(3).
% Output:  The L2 geodesic mean is the rotation matrix S.  

[x1,x2,n] = size(R);
if x1 ~= 3 || x2 ~= 3
   fprintf('The set of rotation matrices should be in SO(3)\n');
   return;
end

S = R(:,:,1);
r = zeros(1,3);
precision = 10^(-5);

while 1
    for i = 1:n
        r = r + rotationMatrixToVector(S'*R(:,:,i)); 
    end
    r = (1/n)*r;
    if norm(r) < precision
        break;
    end
    S = S*rotationVectorToMatrix(r);
end
