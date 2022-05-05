function multiple_perfomance_evaluation()

p_1 = 0.4;%[0,0.1,0.2,0.3,0.4];         % missing data.
p_2 = [0,0.05,0.1,0.15];%[0,0.1,0.2,0.3,0.4];        % outlier.
p_3 = 0.3; %[0,0.2,0.4,0.6,0.8,1];     %noisy arc

p = p_2;
% disp(p);

x_size = size(p,2);
n_rotations = 20;
trials = 5;
err_mean_L1 = zeros(1,x_size);
err_mean_L2 = zeros(1,x_size);

for n=1:x_size
    for trial=1:trials
        fprintf('%d %d\n',n,trial);
        for i=1:n_rotations
            R(3*i-2:3*i,:) = randRotationMatrix();
        end
        M = generating_M_matrix(R,p_1,p_2(n),p_3);
        W_L1 = AvMultipleRotations(1,M);
        W_L2 = AvMultipleRotations(2,M);
        [~,e_L1] = FindRealSolution(R,W_L1,n_rotations);
        [~,e_L2] = FindRealSolution(R,W_L2,n_rotations);
        err_mean_L1(n) = err_mean_L1(n) + e_L1;
        err_mean_L2(n) = err_mean_L2(n) + e_L2;
        
        clear M W_L1 W_L2 e_L1 e_L2 R
    end      
    err_mean_L1(n) = err_mean_L1(n)/trials;
    err_mean_L2(n) = err_mean_L2(n)/trials;
end

figure
plot(p,err_mean_L1,'b--o',p,err_mean_L2,'c--*')
xlabel('Probability for outliers')
ylabel('Error (frobenius norm)')
legend({'L1 AvMultipleRotations','L2 AvMultipleRotations'},'Location','northeast')



