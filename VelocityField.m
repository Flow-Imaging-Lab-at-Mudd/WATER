classdef VelocityField < handle
    properties
        % Imported fields.
        
        % 4D matrix of position vectors on grid.
        X
        % 4D matrix of velocity vectors on grid.
        U
        % Grid dimension.
        dim
        % Optional 4D matrix of noise in velocity to be superposed with U.
        N
        
        % Derived fields.
        
        % Lower and upper bounds of x values.
        xbounds
        % Resolution of x.
        xresol
        % Lower and upper bounds of y values.
        ybounds
        % Resolution of y.
        yresol
        % Lower and upper bounds of z values.
        zbounds
        % Resolution of z.
        zresol
        
    end
    
    methods(Static)
        % Flavors of import functions for different PIV data format.
        
        function VF = import_grid_separate(xw, yw, zw, uw, vw, ww)
           
            % Pack positions and velocities compactly as 3-vectors in extra dimension.
            X = xw;
            X(:,:,:,2) = yw;
            X(:,:,:,3) = zw;
            U = uw;
            U(:,:,:,2) = vw;
            U(:,:,:,3) = ww;
            
            VF = VelocityField(X, U);
        end
        
    end
    
    methods
        % Constructor
        function obj = VelocityField(X, U)
            obj.X = X;
            obj.U = U;
            if ~isequal(size(X), size(U))
                error('Mismatching grid dimensions for position and velocity data')
            end
            dim = size(X);
            obj.dim = dim(1:3);
            % Init null additional noise.
            N = zeros(size(U));
            obj.xbounds = [X(1,1,1,1) X(1,end,1,1)];
            % Suppose our data is not planar.
            obj.xresol = X(1,2,1,1) - X(1,1,1,1);
            obj.ybounds = [X(1,1,1,2) X(end,1,1,2)];
            obj.yresol = X(2,1,1,2) - X(1,1,1,2);
            obj.zbounds = [X(1,1,1,3) X(1,1,end,3)];
            obj.zresol = X(1,1,2,3) - X(1,1,1,3);
        end
        
        % Introduce species of noise.
        function N = noise_uniform(vf, mag, in_place)
            N = rand(size(vf.U))*mag/sqrt(3);
            if exist('in_place', 'var')
                vf.N = N;
                vf.plotNoisyVelocity();
            end
        end
        
        % Plotters
        function plt = plotNoisyVelocity(vf, range)
            if ~exist('range', 'var')
                range = [ones(3, 1) vf.dim'];
            end
            plt = plotVF(vf.X, vf.U + vf.N, range);
        end
        
        % Plotters
        function plt = plotVelocity(vf, range)
            if ~exist('range', 'var')
                range = [ones(3, 1) vf.dim'];
            end
            plt = plotVF(vf.X, vf.U, range);
        end
        
        
    end
end