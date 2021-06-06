classdef VelocityField < handle
    properties
        % Imported fields.
        
        % 4D matrix of position vectors on grid.
        X
        % 4D matrix of velocity vectors on grid.
        U
        % Grid dimension. Note that the indices on grid, due to meshgrid
        % convension, is given as (y, x, z), where x represent the number
        % of distinct values of x, etc..
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
        
        % y and z value may progres from positive to negative.
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
        
        function vf = import_grid_separate(xw, yw, zw, uw, vw, ww)
           
            % Pack positions and velocities compactly as 3-vectors in extra dimension.
            X = xw;
            X(:,:,:,2) = yw;
            X(:,:,:,3) = zw;
            U = uw;
            U(:,:,:,2) = vw;
            U(:,:,:,3) = ww;
            
            vf = VelocityField(X, U);
        end
        
        % Helper for vectorized indices on 4D array.
        function v = getVector(V, index)
            v = squeeze(V(index(2), index(1), index(3), :));
        end
        
    end
    
    methods
        % Constructor taking in valid 4D matrices of position and velocity.
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
            % Suppose our data is not planar and uniformly spaced in
            % position.
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
        
        function ind = getIndices(vf, pos)
            ind = [vf.getIndex_x(pos(1)), vf.getIndex_y(pos(2)), vf.getIndex_z(pos(3))];
        end
        
        
        % Introduce species of noise. Noises generated are added to the
        % present level of noise, not replacing it.
        
        function N = noise_uniform(vf, mag, range)
            if ~exist('range', 'var')
                range = [ones(3, 1) vf.dim'];
            end
            N = zeros(size(vf.U));
            dims = (range(:,2) - range(:,1) + 1)';
            N(range(1,1):range(1,2), range(2,1):range(2,2), range(3,1):range(3,2), :) = ...
                rand([dims 3])*mag/sqrt(3);
            vf.N = vf.N + N;
            vf.plotVelocity(1);
        end
        
        function N = noise_wgn(vf, sd, snr, range)
            if ~exist('range', 'var')
                range = [ones(3, 1) vf.dim'];
            end
            if sd == 0
                N = vf.N;
                vf.N(range(1,1):range(1,2), range(2,1):range(2,2), range(3,1):range(3,2), :) = ...
                    awgn(vf.N(range(1,1):range(1,2), range(2,1):range(2,2), range(3,1):range(3,2), :), snr);
                N = vf.N - N;
            else
                N = zeros(size(vf.U));
                dims = (range(:,2) - range(:,1) + 1)';
                N(range(1,1):range(1,2), range(2,1):range(2,2), range(3,1):range(3,2), :) = ...
                    sd/sqrt(3)*randn([dims 3]);
                vf.N = vf.N + N;
            end
            vf.plotVelocity(1);
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
            title('Velocity $\vec{u}$', 'Interpreter', 'latex')
        end
        
        % The two plotPlane functions can be used to plot any vector field
        % over the grid.
        function plt = plotPlaneSkewed(vf, V, x, eq, with_noise, title_str, range)
            % Plots an arbitrary plane in 3D space given either three
            % non-colinear points or a normal vector + base position paris
            % formula.
            
            if ~exist('range', 'var')
                range = [ones(3, 1) vf.dim'];
            end
            % Obtain a matching 4D boolean matrix indicating membership of
            % points on the plane.
            if sum(size(x), 'all') == 0
                onPlane = skewPlaneMatrix(vf.X, eq(:,1), eq(:,2));
            else
                [~, ~, onPlane] = getPlaneEq(x);
                onPlane = onPlane(vf.X);
            end
            
            if with_noise
                plt = plotVF(vf.X, (V + vf.N) .* onPlane, vf.quiverScale, range);
            else
                plt = plotVF(vf.X, V .* onPlane, vf.quiverScale, range);
            end
            title(title_str, 'Interpreter', 'latex')
        end
        
        function plt = plotPlane(vf, V, with_noise, index, title_str, range)
            % index = [0 0 k], where the nonzero index can be
            % at any dimension, whose values specifies the index of the
            % plane. This is in the usual (x, y, z) orientation.
            
            if ~exist('range', 'var')
                range = [ones(3, 1) vf.dim'];
            end
            % Derive equation for plane.
            eq = zeros(3, 2);
            eq(:, 1) = (index/norm(index))';
            % Obtains a position on the plane.
            eq(:, 2) = VelocityField.getVector(vf.X, index + (index==0));
            % index + (index==0) pads one to the zero-valued components to
            % obtain a valid position index.
            
            plt = vf.plotPlaneSkewed(V, [], eq, with_noise, title_str, range);
        end
        
    end
end