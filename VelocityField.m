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
        
        % Adjustable scale of vector length in plot.
        quiverScale = 1;
        
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
        function vf = VelocityField(X, U)
            vf.X = X;
            vf.U = U;
            if ~isequal(size(X), size(U))
                error('Mismatching grid dimensions for position and velocity data')
            end
            dim = size(X);
            vf.dim = dim(1:3);
            % Init null additional noise.
            vf.N = zeros(size(U));
            vf.xbounds = [X(1,1,1,1) X(1,end,1,1)];
            % Suppose our data is not planar.
            vf.xresol = X(1,2,1,1) - X(1,1,1,1);
            vf.ybounds = [X(1,1,1,2) X(end,1,1,2)];
            vf.yresol = X(2,1,1,2) - X(1,1,1,2);
            vf.zbounds = [X(1,1,1,3) X(1,1,end,3)];
            vf.zresol = X(1,1,2,3) - X(1,1,1,3);
        end

        % Convert spatial coordinates to indices. The position is rounded
        % to the nearest coordinate corresponding to an index.
        function i = getIndex_x(vf, x)
            i = round((x - vf.xbounds(1))/vf.xresol) + 1;
        end
        
        function j = getIndex_y(vf, y)
            j = round((y - vf.ybounds(1))/vf.yresol) + 1;
        end
        
        function k = getIndex_z(vf, z)
            k = round((z - vf.zbounds(1))/vf.zresol) + 1;
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
        function plt = plotVelocity(vf, with_noise, range)
            if ~exist('range', 'var')
                range = [ones(3, 1) vf.dim'];
            end
            if with_noise
                plt = plotVF(vf.X, vf.U + vf.N, vf.quiverScale, range);
            else
                plt = plotVF(vf.X, vf.U, vf.quiverScale, range);
            end
        end
        
        function plt = plotVelocityPlane(vf, x, eq, with_noise, range)
            if ~exist('range', 'var')
                range = [ones(3, 1) vf.dim'];
            end
            % Use given normal vector and base position.
            if sum(size(x), 'all') == 0
                onPlane = skewPlaneMatrix(vf.X, eq(0), eq(1));
            else
                [~, ~, onPlane] = getPlaneEq(x);
                onPlane = onPlane(vf.X);
            end
            
            if with_noise
                plt = plotVF(vf.X, (vf.U + vf.N) .* onPlane, vf.quiverScale, range);
            else
                plt = plotVF(vf.X, vf.U .* onPlane, vf.quiverScale, range);
            end
        end
        
    end
end