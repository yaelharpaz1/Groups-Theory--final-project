function [S,e] = FindRealSolution(R,W,N)
%finding the global rotation matrix.

%Input: R- 3Nx3 arbitrary selction for N rotation matrices. 
       %N- number of rotation matrices.
       %W- 3x3xN, the L1 averaging multiple rotation matrices. 
%Output: S- 3x3xN rotation matrices after minimizing the error. 

%evaluation error function: 
R = reshape(R.',3,3,[]);
R = permute(R,[2,1,3]);

A = zeros(3);
for i=1:N
    A = A + (R(:,:,i)'*W(:,:,i));
end
[u,~,v] = svd(A);
Q = v*u.';                      %Q that minimize the error
S = zeros(3,3,N);

%solution:
for i=1:N
    S(:,:,i) = W(:,:,i)*Q;
end

%error calculation with frobenius norm:
e = 0;
for i=1:N
    e = max(e,norm(R(:,:,i)-S(:,:,i),'fro'));
end

