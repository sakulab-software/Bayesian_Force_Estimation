classdef TFM_computation
    %TFM_CONDITION このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties(Access=public)
        youngs_modulus = 1.0
        poisson_ratio = 0.5
    end
    
    properties(SetAccess=private)
        X, Y         % meshgrid points defined over the substrate: (n, m) shaped
        Fx, Fy       % x- and y-components of cell force, respectively: (n, m) shaped
        IBLx, IBLy   % x- and y-components of initial bead locations, respectively: (B, 1) shaped
        BDx, BDy     % x- and y-components of bead displacements, respectively: (B, 1) shaped
        N            % number of meshgrid points === n times m
        B            % number of beads
        G=0          % elastic matrix: (2B, 2N) shaped
    end
    
    methods
        function obj = TFM_computation(X, Y, Fx, Fy, IBLx, IBLy)
            %TFM_CONDITION
            % X, Y: (n, m) shaped; meshgrid points defined over the substrate
            % Fx, Fy: (n, m) shaped; x- and y-components of cell force, respectively
            % IBLx, IBLy: (B, 1) shaped; x- and y-components of initial bead locations, respectively
            if nargin ~= 6
                error('Just six arguments are expected');
            end
            obj.X = X; obj.Y = Y;
            obj.Fx = Fx; obj.Fy = Fy;
            obj.IBLx = IBLx; obj.IBLy = IBLy;
            obj.N = numel(X);
            obj.B = numel(IBLx);
        end
        
        function this = simulate(this)
            %METHOD1 このメソッドの概要をここに記述
            %   詳細説明をここに記述
            if this.G == 0
                this.G = elastic_tensor(this);
            end
            
            BD = this.G * reshape_F(this);
            [this.BDx, this.BDy] = destruct_BD(this, BD);
            
        end
        
        function [disturbed_BDx, disturbed_BDy] = observe(this, noise_radius, seed)
            rng(seed);
            gm = gmdistribution([0, 0], noise_radius * eye(2));
            noise = random(gm, this.B);
            rng('shuffle');
            
%             BD = sqrt(this.BDx .^ 2 + this.BDy .^ 2);
%             disturbed_BDx = this.BDx + BD .* noise(:, 1);
%             disturbed_BDy = this.BDy + BD .* noise(:, 2);
            disturbed_BDx = this.BDx + noise(:, 1);
            disturbed_BDy = this.BDy + noise(:, 2);
        end
        
        function [BDx, BDy] = observe_without_noise(this)
            BDx = this.BDx;
            BDy = this.BDy;
        end
    end
end


function G = elastic_matrix(obj)

G = zeros(2 * obj.B, 2 * obj.N);
nu = obj.poisson_ratio;
coef = (1 + nu) / pi / obj.youngs_modulus;

for b=1:obj.B
    for n=1:obj.N
        x = obj.IBLx(b) - obj.X(n);
        y = obj.IBLy(b) - obj.Y(n);
        r = sqrt(x ^ 2 + y ^ 2);
        tmp = coef / r ^ 3 * ... 
            [(1 - nu) * r ^ 2 + nu * x ^ 2, nu * x * y;
             nu * x * y,                    (1 - nu) * r ^ 2 + nu * y ^ 2];
        G(2 * b - 1:2 * b, 2 * n - 1: 2 * n) = tmp;
    end
end

end

function G = elastic_tensor(obj)

G = zeros(2 * obj.B, 2 * obj.N);
nu = obj.poisson_ratio;
coef = (1 + nu) / pi / obj.youngs_modulus;

x = obj.IBLx - reshape(obj.X, 1, []);
y = obj.IBLy - reshape(obj.Y, 1, []);
r = sqrt(x .^ 2 + y .^ 2);

G11 = coef ./ r .^ 3 .* ((1 - nu) .* r .^ 2 + nu .* x .^ 2);
G12 = coef ./ r .^ 3 .* (nu .* x .* y);
G22 = coef ./ r .^ 3 .* ((1 - nu) .* r .^ 2 + nu .* y .^ 2);

G(1:2:size(G, 1), 1:2:size(G, 2)) = G11;
G(2:2:size(G, 1), 1:2:size(G, 2)) = G12;
G(1:2:size(G, 1), 2:2:size(G, 2)) = G12;
G(2:2:size(G, 1), 2:2:size(G, 2)) = G22;

end

function F = reshape_F(obj)

F = zeros(2 * obj.N, 1);

for n=1:obj.N
    F(2 * n - 1: 2 * n) = [obj.Fx(n); obj.Fy(n)];
end

end

function [BDx, BDy] = destruct_BD(obj, BD)

BDx = zeros(obj.B, 1);
BDy = zeros(obj.B, 1);

for b=1:obj.B
    BDx(b) = BD(2 * b - 1);
    BDy(b) = BD(2 * b);
end

end