function performance_evaluation_noise()

noise_level = [0,1,2,3,4,5,6,7,8,9,10];
n_rotations = 100;
level = size(noise_level,2);
trials = 200;                            
err_mean_L1 = zeros(1,level);
err_mean_L2 = zeros(1,level);

for n=1:level
    for trial=1:trials
        R = randRotationMatrix();           
        axang = rotm2axang(R);
        r = axang(1:3);
        angle = axang(4);
        RR_err = zeros(3,3,n_rotations);
        for i=1:n_rotations
            theta = normrnd(0,n);               % (degrees)
            while theta > n+3                    % No outliers
                theta = normrnd(0,n);           % (degrees)
            end
            theta = theta*pi/180;               % (radians)
            RR_err(:,:,i) = R*axang2rotm([r(1) r(2) r(3) theta]);
        end
        R_L1 = L1_geodesic_mean(RR_err);
        R_L2 = L2_geodesic_mean(RR_err);
        
        axang_L1 = rotm2axang(R_L1);
        angle_L1 = abs(axang_L1(4) - angle);
        err_mean_L1(n) = err_mean_L1(n) + angle_L1;
        
        axang_L2 = rotm2axang(R_L2);
        angle_L2 = abs(axang_L2(4) - angle);
        err_mean_L2(n) = err_mean_L2(n) + angle_L2;
    end
    err_mean_L1(n) = err_mean_L1(n)/trials;
    err_mean_L2(n) = err_mean_L2(n)/trials;
end

% Convert to degrees
err_mean_L1 = err_mean_L1*180/pi;
err_mean_L2 = err_mean_L2*180/pi;

figure
plot(noise_level,err_mean_L1,'b--o',noise_level,err_mean_L2,'c--*')
xlabel('Noise (degrees)')
ylabel('Error (degrees)')
legend({'L1 rotation averaging','L2 rotation averaging'},'Location','northwest')

% 
% % A normally distributed angle θ with mean 0 and standard deviation σ = 2 degrees 
% % is generated to simulate the effect of random rotation noise.
% theta = normrnd(0,2);       % (degrees)
% theta = theta*2*pi/360;
% 
% % no outliers: theta < 5 degrees.
