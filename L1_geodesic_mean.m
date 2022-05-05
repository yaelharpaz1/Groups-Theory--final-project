function S = L1_geodesic_mean(R)

% Computing the L1 geodesic mean in a group of rotation matrices 
% using the Weiszfeld algorithm.

% Input:   R, a 3x3xn matrix, represents a set of rotation matrices in SO(3).
% Output:  The L1 geodesic mean is the rotation matrix S. 

[x1,x2,n] = size(R);
if x1 ~= 3 || x2 ~= 3
   fprintf('The set of rotation matrices should be in SO(3)\n');
   return;
end

v         = zeros(3,n);
norms     = zeros(1,n);
precision = 10^(-10);

% Initial estimation:
S = L2_chordal_mean(R);

while isreal(S)
    
    S_prev = S;
    
    % Apply the logarithm map, Using the Rodrigues formula:
    for i=1:n
        v(:,i)   = rotationMatrixToVector(R(:,:,i)*S^(-1));
        if norm(v(:,i)) ~= 0
            norms(i) = 1/norm(v(:,i));
        else
            norms(i) = 0;
        end
        v(:,i)   = v(:,i)*norms(i);
    end
    
    % Weiszfeld step
    delta = sum(v,2)/sum(norms);
    
    % Apply the exp map, Using the Rodrigues formula:
    S = rotationVectorToMatrix(delta)*S;
    
    % Convergence criterion:
    % The distance between the current rotation matrix and the previous one.
    % The chordal metric:
    if norm(S - S_prev,'fro') < precision
        break;
    end
end



