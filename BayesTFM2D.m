classdef BayesTFM2D
    %BAYESTFM このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties(Access=private)
        alpha, beta, thr % regularization parameters
        max_iter    % number of max iterations
    end
    
    properties(Access=public)
        force
        mse
    end
    
    methods
        function obj = BayesTFM2D(alpha, beta, thr, max_iter)
            if nargin ~= 4
                error('Just four arguments are required');
            end
            
            obj.alpha = alpha;
            obj.beta = beta;
            obj.thr = thr;
            obj.max_iter = max_iter;
        end
        
        function obj = fit(obj, BD, G, B, X, Y, IBLx, IBLy)
            n = length(B) / 2;
            M = magnitude_array(B, 2);
            D = reshape(repmat(M', 2, 1)', 1, [])';
            D = B ./ D;

%             [obj.mse, obj.force] = bayes_estimation(...
%                 n, BD', B', M', D', G, obj.alpha, obj.beta, obj.thr, obj.max_iter);

            [obj.mse, obj.force] = bayes_estimation_with_sampling(...
                n, BD', B', M', D', G, obj.alpha, obj.beta, obj.thr, obj.max_iter, 20, X, Y, IBLx, IBLy);
        end
    end
end

function M = magnitude_array(v, d)
M = zeros(length(v) / d, 1);
for i=1:length(v) / d
    M(i) = norm(v(d * (i - 1) + 1:d * i));
end
end

function [err, estimated_f] = bayes_estimation(...
    division_number, RBD, IE, mag_IE, dir_IE, G, alpha, beta, thr, max_iter)

Scale = eye(length(IE));
Rotation_array = zeros([2, length(IE)]);
for i=1:length(IE) / 2
    Rotation_array(1, 2 * i - 1) = 1;
    Rotation_array(1, 2 * i) = 0;
    Rotation_array(2, 2 * i - 1) = 0;
    Rotation_array(2, 2 * i) = 1;
end

[err, estimated_f] = bayes_estimation_(...
    division_number, RBD, IE, mag_IE, dir_IE, G, alpha, beta, thr, max_iter, Scale, Rotation_array);

end

function [err, estimated_f] = bayes_estimation_with_sampling(...
    division_number, RBD, IE, mag_IE, dir_IE, G, alpha, beta, thr, max_iter, samplingNumber, X, Y, IBLx, IBLy)

Cerr = inf(1, samplingNumber);
Cest = cell(1, samplingNumber);

for i=1:samplingNumber
    rng(i + 1);
%     rng("shuffle");

%     Scale = gamrnd(2, 1, [1, length(mag_IE)]);
%     Scale = dimeye(Scale, 2);
    Scale = eye(length(IE));

    [nBDx, nBDy] = nearestBeadDisp(X, Y, IBLx, IBLy, RBD(1:2:end), RBD(2:2:end));
    Rotation_array = rotateMatrix(dir_IE(1:2:end), dir_IE(2:2:end), nBDx, nBDy, 1e-1 * randn(1, division_number ^ 2));

%     Rotation_array = zeros([2, length(IE)]);
%     for j=1:length(IE) / 2
%         Rotation_array(1, 2 * i - 1) = 1;
%         Rotation_array(1, 2 * i) = 0;
%         Rotation_array(2, 2 * i - 1) = 0;
%         Rotation_array(2, 2 * i) = 1;
%     end
    
    [Cerr(i), Cest{i}] = bayes_estimation_(...
        division_number, RBD, IE, mag_IE, dir_IE, G, alpha, beta, thr, max_iter, Scale, Rotation_array);
end

[err, I] = min(Cerr);
estimated_f = Cest{I};

end

function [err, estimated_f] = bayes_estimation_(...
    division_number, RBD, IE, mag_IE, dir_IE, G, alpha, beta, thr, max_iter, initialS, initialR)

Scale = initialS;
Rotation_array = initialR;

Sigma_f_inv = alpha * (G' * G) + beta * eye(length(IE));
% Sigma_f = inv(Sigma_f_inv);
alphaSigma_fGU = alpha * RBD * G / Sigma_f_inv;

err = mean((RBD - IE * G') .^ 2);
estimated_f = IE;

prev_err = inf;
prev_estimated_f = ones([1 2 * division_number ^ 2]);

mu_f = zeros([1, length(IE)]);
mag_mu_f = zeros([1, length(IE) / 2]);

for i=1:max_iter
    if i > max_iter
        disp("may not converged");
    end
    % Estep
    mu_f = alphaSigma_fGU + beta * IE * Scale * block_diag(Rotation_array)' / Sigma_f_inv;
    
    mu_f(mu_f > thr) = mu_f(mu_f > thr);% - thr;
    mu_f(mu_f < -thr) = mu_f(mu_f < -thr);% + thr;
    mu_f(mu_f <= thr & mu_f >= -thr) = ...
        0 * sign(mu_f(mu_f <= thr & mu_f >= -thr)); % modified: 1e12 -> 0
    
    % Mstep
    for j=1:length(IE) / 2
        mag_mu_f(j) = norm([mu_f(2 * j - 1), mu_f(2 * j)]);
%         mag_mu_f(j) = 1 / 100 * norm([mu_f(2 * j - 1), mu_f(2 * j)]);
    end
        
    for j=1:length(IE) / 2
       if mag_IE(j) > 0
           theta1 = angle(complex(mu_f(2 * j - 1), mu_f(2 * j)));
           theta2 = angle(complex(dir_IE(2 * j - 1), dir_IE(2 * j)));
           theta = theta1 - theta2;
           cosine_value = cos(theta);
           sine_value = sin(theta);
           Rotation_array(1, 2 * j - 1) = cosine_value;
           Rotation_array(1, 2 * j) = - sine_value;
           Rotation_array(2, 2 * j - 1) = sine_value;
           Rotation_array(2, 2 * j) = cosine_value;
           
           Scale(2 * j - 1, 2 * j - 1) = mag_mu_f(j) / mag_IE(j);
           Scale(2 * j, 2 * j) = Scale(2 * j - 1, 2 * j - 1);
       end
    end
    prev_err = err;
    prev_estimated_f = estimated_f;
    estimated_f = IE * Scale * block_diag(Rotation_array)';
    err = mean((RBD - estimated_f * G') .^ 2) * 2;
    
    if prev_err < err
%         err = prev_err;
%         estimated_f = prev_estimated_f;
%         break;
    elseif (prev_err - err) / prev_err < 10 ^ -3
        break;
    end
end

end

function [M] = dimeye(v, k)

v = reshape(v, 1, []);
M = diag(reshape(repmat(v, [k, 1]), [], 1)');

end

function [M] = block_diag(array)

[dimension_of_one_block, number_of_blocks] = size(array);
number_of_blocks = number_of_blocks / dimension_of_one_block;

M = zeros(number_of_blocks * dimension_of_one_block);
for i=0:number_of_blocks - 1
    for j=1:dimension_of_one_block
        for k=1:dimension_of_one_block
            M(dimension_of_one_block * i + j, dimension_of_one_block * i + k) = array(j, dimension_of_one_block * i + k);
        end
    end
end
end

function [nBDx, nBDy] = nearestBeadDisp(X, Y, BX, BY, BDx, BDy)

N = numel(X);
B = numel(BX);

Dist = sqrt((X(:) - BX(:)') .^ 2 + (Y(:) - BY(:)') .^ 2);
[~, I] = min(Dist, [], 2);

nBDx = zeros(size(X));
nBDy = zeros(size(X));

for i=1:N
    m = sqrt(BDx(I(i)) ^ 2 + BDy(I(i)) ^ 2);
    nBDx(i) = BDx(I(i)) / m;
    nBDy(i) = BDy(I(i)) / m;
end

end

function [R] = rotateMatrix(Ax, Ay, Bx, By, r)

R = zeros(2, 2 * numel(Ax));
for i=1:numel(Ax)
   theta1 = angle(complex(Bx(i), By(i)));
   theta2 = angle(complex(Ax(i), Ay(i)));
   theta = theta1 - theta2 + r(i);
   cosine_value = cos(theta);
   sine_value = sin(theta);
   R(1, 2 * i - 1) = cosine_value;
   R(1, 2 * i) = - sine_value;
   R(2, 2 * i - 1) = sine_value;
   R(2, 2 * i) = cosine_value;
end

end
