function M = generating_M_matrix(R,p_1,p_2,p_3)

%generating M matrix for demonstrate purpose.
%noise/missing data analysis- generalization over SO(3).

%Input: R- 3Nx3 arbitrary selction for N rotation matrices. 
%Output: M- 3Nx3N block matrix of relative rotation matrices Rij.

m = size(R,1); 
N = m/3;                %number of unknown rotation matrices.

U = R*R.';              %3Nx3N block matrix of relative rotation matrices Rij without noise.

%creating M matrix:
M = diag(ones(3*N,1)); 

for i=1:N
    for j=(i+1):N
        x = rand;
        if x < p_1
            M(3*i-2:3*i,3*j-2:3*j) = zeros(3);          %missing data.
        else
            x = rand;
            if x < p_2
                M(3*i-2:3*i,3*j-2:3*j) = randRotationMatrix();      %random matrix in SO(3).
            else
                x = rand;
                if x < p_3
                    J = U(3*j-2:3*j,3*i-2:3*i);
                    axang = rotm2axang(J);
                    r = axang(1:3);
                    theta = normrnd(0,2);               % (degrees)
                    while theta > 4                     % No outliers
                        theta = normrnd(0,2);           % (degrees)
                    end
                    theta = theta*pi/180;               % (radians)
                    J = J*axang2rotm([r(1) r(2) r(3) theta]);
                    M(3*i-2:3*i,3*j-2:3*j) = J;
                else
                    M(3*i-2:3*i,3*j-2:3*j) = U(3*j-2:3*j,3*i-2:3*i);
                end
            end 
        end
    end
end

M = triu(M) + tril(M.',-1); %enforce symmetry.

