function R = AvMultipleRotations(Lp,M)

%minimizing the cost function of R using geodesic distances.
%using the Weiszfeld algorithm.
%assumption: for a large data set, each i has at least one neighbour j.
%excellent results on large data set.

%Input: M- 3Nx3N block matrix of relative rotation matrices Rij, that provided
        %by some measurement process (to demonstrate use generating_M_matrix algorithm)        
%Output: R- 3x3xN, the L1 averaging multiple rotation matrices. 

m = size(M,1);
N = m/3;                          % number of unknown rotation matrices.
R = zeros(3,3,N);                 % the rotation matrices.
precision = 10^(-10);

%1. Initialization (set Ri0 with the maximum number of neighbours to the identity rotation):

x = 0;
for i=1:N
    z = nnz(M(i+2*(i-1):3*i,:));
    if z > x
        x  = z;
        i0 = i;
    end
end
R(:,:,i0) = eye(3);

R_es = zeros(3*N,3,N);              %estimations for each Ri

%Ri (i~=i0) estimation from Ri0 (by Ri=Ri0i*Ri0):
for i=1:N
    if i ~= i0
        if nnz(M(3*i0-2:3*i0,3*i-2:3*i)) ~= 0
            R_es(3*i0-2:3*i0,:,i) = M(3*i0-2:3*i0,3*i-2:3*i);
        end
    end
end

%Ri (i~=i0) estimations from Rj (j~=io)(by Ri=Rji*Rj):
for i=1:N
    if i ~= i0
        for j=1:N
            if j ~= i0 && j ~= i
                if nnz(M(3*j-2:3*j,3*i-2:3*i)) ~= 0
                    R_es(3*j-2:3*j,:,i) = M(3*j-2:3*j,3*i-2:3*i)*R_es(3*i0-2:3*i0,:,j);
                end
            end
        end
    end
end

%2. sweep:
while 1
    R_prev = R;
    %for each i- computing the L1 geodesic mean by one iteration of the Weiszfeld algorithm: 
    for i=1:N
        if i ~= i0  
            T = R_es(:,:,i);
            T = T(any(T,2),:);                          %remove zeros rows 
            T = reshape(T.',3,3,[]);
            T = permute(T,[2,1,3]);
            if Lp == 1
                R(:,:,i) = L1_geom_1iteration_of_WA(T);
            end
            if Lp == 2
                R(:,:,i) = L2_chordal_mean(T);
            end
        end
    end
    
    %convergence criterion:
    d = 0;
    for t=1:N
        d = max(norm(R(:,:,t) - R_prev(:,:,t),'fro'),d);
    end
    if d < precision
        break;
    end
    
    %Ri (i~=i0) estimations from Rj (by Ri=Rji*Rj):
    R_es = zeros(3*N,3,N);
    for i=1:N
        if i ~= i0
            for j=1:N
                if j ~= i
                    if nnz(M(3*j-2:3*j,3*i-2:3*i)) ~= 0
                        R_es(3*j-2:3*j,:,i) = M(3*j-2:3*j,3*i-2:3*i)*R(:,:,j);
                    end
                end
            end
        end
    end
end


        
        
  


       
            
        





