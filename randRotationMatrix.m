function M = randRotationMatrix()

x1 = 2*pi*rand;
x2 = 2*pi*rand;
x3 = rand;
R = [cos(x1),sin(x1),0;-sin(x1),cos(x1),0;0,0,1];
v = [cos(x2)*sqrt(x3),sin(x2)*sqrt(x3),sqrt(1-x3)]';
H = eye(3) - 2*v*v';                                    % Householder matrix
M = -H*R;                                               % M is uniformly distributed within SO(3)    
