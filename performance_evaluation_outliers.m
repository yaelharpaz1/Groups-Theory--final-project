function performance_evaluation_outliers()

outliers_pre = [0,5,10,15,20,25,30,35,40];
x_size = size(outliers_pre,2);
n_rotations = 100;
trials = 200;                            
err_mean_L1 = zeros(1,x_size);
err_mean_L2 = zeros(1,x_size);

for n=1:x_size
    for trial=1:trials
        R = randRotationMatrix();           
        axang = rotm2axang(R);
        r = axang(1:3);
        angle = axang(4);
        RR_err = zeros(3,3,n_rotations);
        for i=1:outliers_pre(n)
            theta = normrnd(0,20);
            while theta < 5                     
                theta = normrnd(0,20);           % (degrees)
            end
            theta = theta*pi/180;               % (radians)
            RR_err(:,:,i) = R*axang2rotm([r(1) r(2) r(3) theta]);
        end
        for i=(outliers_pre(n)+1):n_rotations
            theta = normrnd(0,2);               % (degrees)
            while theta > 5                     % No outliers
                theta = normrnd(0,2);           % (degrees)
            end
            theta = theta*pi/180;               % (radians)
            RR_err(:,:,i) = R*axang2rotm([r(1) r(2) r(3) theta]);
        end
        R_L1 = L1_geodesic_mean(RR_err);
        R_L2 = L2_chordal_mean(RR_err);
        
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
plot(outliers_pre,err_mean_L1,'b--o',outliers_pre,err_mean_L2,'c--*')
xlabel('Outlier percentage (%)')
ylabel('Error (Degrees)')
legend({'L1 rotation averaging','L2 rotation averaging'},'Location','northwest')
