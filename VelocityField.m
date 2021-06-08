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
        dims
        % Optional 4D matrix of noise in velocity to be superposed with U.
        N
        
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
        % Array of resolutions.
        resol
        
        % Range for computation and visualization.
        range
        % Length of range in 3 dimensions.
        span
        
        % struct for properties of solvers.
        solver
        % struct for graphers.
        plotter
        % struct for invariant attributes of fluid.
        fluid
        
    end
    
    methods(Static)
        
        %%%%%%%%%%%%%%%%%%% Adapters of Recorded Data %%%%%%%%%%%%%%%%%%%
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
        
        %%%%%%%%%%%%%%%%%%% Various Helpers %%%%%%%%%%%%%%%%%%%%%%
        function v = getVector(V, index)
            % Helper for vectorized indices on 4D array.
            % Non-vectorized; index can only be a vector.
            v = squeeze(V(index(2), index(1), index(3), :));
        end
        
        function v = getV(V, index)
        end
        
        % Helpers for computing errors given reference.
        
        function err = error_L2(V, V0)
            err = sum((V - V0).^2, 4);
        end
        
    end
    
    methods
        
        function vf = VelocityField(X, U)
            % Constructor taking in valid 4D matrices of position and
            % velocity whose first three dimensions conform to that of a
            % meshgrid.
            
            vf.X = X;
            vf.U = U;
            if ~isequal(size(X), size(U))
                error('Mismatching grid dimensions for position and velocity data')
            end
            dims = size(X);
            vf.dims = dims(1:3);
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
            vf.resol = [vf.xresol vf.yresol vf.zresol];
            vf.range = [ones(3, 1) vf.dims'];
            vf.span = (vf.range*[-1; 1])' + 1;
            
            vf.initPropertyStructs();
            set(0,'defaultTextInterpreter','latex');
        end
        
        function initPropertyStructs(vf)
            % Default innate attributes.
            vf.fluid.density = 1;
            
            % Default solver attributes.
            vf.solver.dv = abs(vf.xresol*vf.yresol*vf.zresol);
            
            % Properties of differentiation.
            vf.solver.diff.order = 1;
            vf.solver.diff.mode = 'central';
            
            vf.solver.ke.mode = 'direct';
            
            % Default plotter attributes.
            vf.plotter.quiverScale = 1;
        end

        %%%%%%%%%%%%%%%%%%%%% Coordinate Helpers %%%%%%%%%%%%%%%%%%%%%
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
        
        function x = get_x(vf, i)
            x = vf.xbounds(1) + (i-1)*vf.xresol;
        end
        
        function y = get_y(vf, j)
            y = vf.ybounds(1) + (j-1)*vf.yresol;
        end
        
        function z = get_z(vf, k)
            z = vf.zbounds(1) + (k-1)*vf.zresol;
        end
        
        function eq = getRegPlaneEq(vf, index)
            % For a regular plane, one which is perpendicular to one of the
            % basis unit vectors, i.e. x y z, obtain the standard
            % representation as orthogonal vector and base position.
            %
            % For example, index = [0 0 k]', where k is the index of the
            % plane in the appropriate dimension.
            %
            % Non-vectorized, index a 3 x 2 matrix.
            
            eq = zeros(3, 2);
            eq(:, 1) = index/sum(index, 1);
            % Obtain a position on the plane.
            eq(:, 2) = VelocityField.getVector(vf.X, index + (index==0));
        end
        
        function eqs = getRegPlaneEqs(vf, indices)
            % Index vectors stacked as columns in a 2D array.
            
            eqs = zeros(size(indices, 2), 3, 2);
            for i = 1: size(indices, 2)
                eqs(i, :, :) = vf.getRegPlaneEq(indices(:, i));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%% Add Noise %%%%%%%%%%%%%%%%%%%%%%
        
        % Introduce species of noise. Noises generated are added to the
        % present level of noise, not replacing it.
        function N = noise_uniform(vf, mag, range)
            if ~exist('range', 'var')
                range = [ones(3, 1) vf.dims'];
            end
            N = zeros(size(vf.U));
            dims = (range(:,2) - range(:,1) + 1)';
            N(range(1,1):range(1,2), range(2,1):range(2,2), range(3,1):range(3,2), :) = ...
                rand([dims 3])*mag/sqrt(3);
            vf.N = vf.N + N;
            vf.plotVector(vf.U, N, '$\vec{u}$');
        end
        
        function N = noise_wgn(vf, sd, snr, range)
            % Add white gaussian noise specified by either a standard
            % deviation or a signal-to-noise ratio.
            
            if ~exist('range', 'var')
                range = [ones(3, 1) vf.dims'];
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
            vf.plotVector(vf.U, N, '$\vec{u}$');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%% Plotters %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function plt = plotVector(vf, V, noise, title_str, range)
            if ~exist('range', 'var')
                range = [ones(3, 1) vf.dims'];
            end
            
            plt = plotVF(vf.X, V + noise, vf.plotter.quiverScale, range);
            
            title(title_str)
        end
        
        function plt = plotScalar(vf, S, noise, title_str, range)
            % Show a color map of the scalar field over the range
            % specified, defaulted to global.
            
            if ~exist('range', 'var')
                range = [ones(3, 1) vf.dims'];
            end
            % Subset region plotted.
            X = vf.X(range(1,1): range(1,2), range(2,1): range(2,2), range(3,1): range(3,2), :);
            S = S(range(1,1): range(1,2), range(2,1): range(2,2), range(3,1): range(3,2), :);
            % Collapse 4D matrix into a set of points.
            X = reshape(X, [], 3);
            S = S(:);
            % Scaled size of each dot displayed. TODO: more fitting
            % expression.
            dif = sort(range(:,2) - range(:,1), 2);
            dot_size = 40^(sum(dif ~= 0)-1)/(dif(1)*dif(2));
            
            plt = figure;
            scatter3(X(:,1), X(:,2), X(:,3), dot_size, S + noise, 'filled');
            
            colorbar
            xlabel('$x$')
            ylabel('$y$')
            zlabel('$z$')
            title(title_str)
        end
        
        function plt = slicePlanes(vf, S, noise, planes, title_str, range)
            % SLICEPLANES renders continuous color plots on the planes
            % whose equations are specified. planes are in the format of
            % [[orth base]...] for all the planes to be plotted.
            
            plt = figure;
            % Subset grid.
            if ~exist('range', 'var')
                range = [ones(3, 1) vf.dims'];
            end
            X = vf.X(range(1,1): range(1,2), range(2,1): range(2,2), range(3,1): range(3,2), :);
            S = S(range(1,1): range(1,2), range(2,1): range(2,2), range(3,1): range(3,2), :);
            
            for i = 1: size(planes, 1)
                eq = squeeze(planes(i, :, :));
                z = (dot(eq(:,2), eq(:,1)) - eq(1,1)*X(:,:,1,1) - eq(2,1)*X(:,:,1,2)) / eq(3,1);
                slice(squeeze(X(:,:,:,1)), squeeze(X(:,:,:,2)), squeeze(X(:,:,:,3)), S + noise, ...
                    squeeze(X(:,:,1,1)), squeeze(X(:,:,1,2)), z)
                hold on
            end
            
            colorbar
            xlabel('$x$')
            ylabel('$y$')
            zlabel('$z$')
            title(title_str)
        end
        
        
        % The two plotPlane functions can be used to plot any vector field
        % over the grid.
        function plt = plotPlaneVector(vf, V, eq, noise, title_str, range)
            % Plots an arbitrary plane in 3D space given either three
            % non-colinear points or a normal vector + base position paris
            % formula.
            
            if ~exist('range', 'var')
                range = [ones(3, 1) vf.dims'];
            end
            % Obtain a matching 4D boolean matrix indicating membership of
            % points on the plane.
            onPlane = skewPlaneMatrix(vf.X, eq(:,1), eq(:,2), 3);
            
            plt = plotVF(vf.X, (V + noise) .* onPlane, vf.plotter.quiverScale, range);
            title(title_str)
        end
        
        function plt = plotScalarPlaneSkewed(vf, Mag, noise, eq, title_str)
            % Derive plane equation if a normal plane is given.
            if isequal(eq(:, 2), zeros(3, 1))
                
            end
            % Boolean predicate for determining if points are on plane.
            onPlane = skewPlaneMatrix(vf.X, eq(:, 1), eq(:, 2), 1);
        end
        
        function plt = plotPlaneScalar(vf, S, noise, range, title_str)
            % Plots a scalar field Mag over X.
            % Mag is a 3D matrix of corresponding to positions in X.
            % Only a regular plane is currently allowed.
            
            % Find dimension, i.e. x,y,z, perpendicular to the plane.
            index = find((range(:,2) - range(:,1)) == 0, 1);
            % Extract the plane.
            x = squeeze(vf.X(range(1,1): range(1,2), range(2,1): range(2,2), range(3,1): range(3,2), :));
            S = squeeze(S(range(1,1): range(1,2), range(2,1): range(2,2), range(3,1): range(3,2)));
            % Dimensionalize noise to allow scalar input.
            if isequal(size(noise), [1 1])
                noise = repmat(noise, vf.dims);
            end
            noise = squeeze(noise(range(1,1): range(1,2), range(2,1): range(2,2), range(3,1): range(3,2)));
            
            plt = figure;
            
            switch index
                case 1 %y
                    % Non right-handed space.
                    surf(x(:,:,1), x(:,:,3), S + noise);
                    xlabel('$x$')
                    ylabel('$z$')
                case 2 %x
                    surf(x(:,:,2), x(:,:,3), S + noise);
                    xlabel('$y$')
                    ylabel('$z$')
                case 3 %z
                    surf(x(:,:,1), x(:,:,2), S + noise);
                    xlabel('$x$')
                    ylabel('$y$')
            end
            title(title_str)
        end
        
        function plt = isosurfaces(vf, S, isovals, title_str, range)
            % Given an array of values of a scalar field S, their
            % isosurfaces are plotted.
            
            if ~exist('range', 'var')
                range = [ones(3, 1) vf.dims'];
            end
            % Subset range.
            X = squeeze(vf.X(range(1,1): range(1,2), range(2,1): range(2,2), range(3,1): range(3,2), :));
            S = squeeze(S(range(1,1): range(1,2), range(2,1): range(2,2), range(3,1): range(3,2)));
            
            plt = figure;
            for isoval = isovals
                isosurface(X(:,:,:,1), X(:,:,:,2), X(:,:,:,3), S, isoval)
                hold on
            end
            colorbar
            xlabel('$x$')
            ylabel('$y$')
            zlabel('$z$')
            title(title_str)
        end
        
        %%%%%%%%%%%%%%%%%%% Solvers of Derived Quantities %%%%%%%%%%%%%%%%%%
        
        function k = kineticEnergy(vf, range)
            
            % Subset range of interest.
            if ~exist('range', 'var')
                range = [ones(3, 1) vf.dims'];
            end
            U = squeeze(vf.U(range(1,1): range(1,2), range(2,1): range(2,2), range(3,1): range(3,2), :));
            
            % Selected mode of computation.
            switch vf.solver.ke.mode
                case 'direct'
                    k = 1/2*vf.fluid.density*vf.solver.dv*sum(U.^2, 'all');
            end
        end
        
        %%%%%%%%%%%%%%%%%%% Differential Methods %%%%%%%%%%%%%%%%%%%%%%
        
        function nder = diff1(vf, V, ind_inc, mode, range)
            % nder = vf.DIFF1(V, ind_inc, mode) computes the first order
            % (simple difference) numerical direction derivative of 'V' as
            % specified by the vector increment in index 'ind_inc'. The
            % parameter 'mode' specified the finite difference scheme used:
            % 'right', 'left', or 'central'. The derivative is computed in
            % the region of the volume for which the difference exists.

            % Subset for region of interest.
            if ~exist('range', 'var')
                range = [ones(3, 1) vf.dims'];
            end
            V = squeeze(V(range(1,1): range(1,2), range(2,1): range(2,2), range(3,1): range(3,2), :));
            
            % Range of differentiable values, from either corners, unspecified.
            rem = size(V, 1:3) - [ind_inc(2) ind_inc(1) ind_inc(3)];

            switch mode
                case 'right'
                    right_shifted = NaN(size(V));
                    % Shift according to direction and increment in index.
                    right_shifted(1:rem(1), 1:rem(2), 1:rem(3), :) = ...
                        V(ind_inc(2)+1: end, ind_inc(1)+1: end, ind_inc(3)+1: end, :);
                    dif = right_shifted - V;
                case 'left'
                    left_shifted = NaN(size(V));
                    % Shift according to direction and increment in index.
                    left_shifted(ind_inc(2)+1: end, ind_inc(1)+1: end, ind_inc(3)+1: end, :) = ...
                        V(1:rem(1), 1:rem(2), 1:rem(3), :);
                    dif = V - left_shifted;
                case 'central'
                    nder = (vf.diff1(V, ind_inc, 'right', range) + vf.diff1(V, ind_inc, 'left', range)) / 2;
                    return
            end
            % Compute displacement in direction and differentiate.
            nder = dif / norm(vf.resol .* ind_inc);
        end
        
        function jacob = jacobian(vf, V, mode)
            % Computes the jacobian matrix (or gradient, or derivative
            % matrix) of 'V' with numerical scheme given in 'mode'.
            
            jacob = NaN([size(V, 1:3) 3 size(V, 4)]);
            
            % Subset region of interest.
            if ~isequal(size(V, 1:3), span)
                V = squeeze(V(vf.range(1,1): vf.range(1,2), vf.range(2,1): vf.range(2,2), ...
                        vf.range(3,1): vf.range(3,2), :));
            end
            
            switch mode
                case 'unit'
                    jacob(:, :, :, 1, :) = vf.diff1(V, [1 0 0], vf.solver.diff.mode);
                    jacob(:, :, :, 2, :) = vf.diff1(V, [0 1 0], vf.solver.diff.mode);
                    jacob(:, :, :, 3, :) = vf.diff1(V, [0 0 1], vf.solver.diff.mode);
            end
        end
        
        function div = div(vf, V)
            div = NaN(size(V, 1:3));
             % Subset region of interest.
            if ~isequal(size(V, 1:3), vf.span)
                V = squeeze(V(vf.range(1,1): vf.range(1,2), vf.range(2,1): vf.range(2,2), ...
                        vf.range(3,1): vf.range(3,2), :));
            end
            
            switch vf.solver.diff.order
                case 1
                    % Ensure diff direction in increasing x, y, z by boolean multiplication.
                    div = squeeze(vf.diff1(V(:,:,:,1), [1 0 0], vf.solver.diff.mode)*(vf.xresol>=0) + ...
                        vf.diff1(V(:,:,:,2), [0 1 0], vf.solver.diff.mode)*(vf.yresol>=0) + ...
                        vf.diff1(V(:,:,:,3), [0 0 1], vf.solver.diff.mode)*(vf.zresol>=0));
            end
        end
        
    end
end