function S = L2_chordal_mean(R)

% Input:   R, a 3x3xn matrix, represents a set of rotation matrices in SO(3).
% Output:  The L2 chordal mean is the rotation matrix S.  

[x1,x2,n] = size(R);
if x1 ~= 3 || x2 ~= 3
   fprintf('The set of rotation matrices should be in SO(3)\n');
   return;
end

S = zeros(3,3);

for i = 1:n
   S = S + R(:,:,i); 
end

[U,~,V] = svd(S);

if det(U*V') >= 0
    S = U*V';
else
    S = U*diag([1 1 -1])*V';
end