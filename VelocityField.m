% Object representation of a 3D PIV Velocity Field.
%
% Derek Li & Leah Mendelson

classdef VelocityField < handle
    properties
        % Imported fields.
        
        % 4D matrix of position vectors on grid.
        X
        % 4D matrix of velocity vectors on grid.
        U
        
        % 4D matrix of vorticity vectors derived from velocity.
        vort
        
        % X, U subsetted into region of interest by vf.range.
        X_e
        U_e
        vort_e
        
        % Grid dimension. Note that the indices on grid, due to meshgrid
        % convension, is given as (y, x, z), where x represent the number
        % of distinct values of x, etc.. Row vector.
        dims
        % Optional 4D matrix of noise in velocity to be superposed with U.
        N
        % Noise in the region of interest.
        N_e
%         % Vortivity noise.
%         Nvort
%         Nvort_e
        
        % Commonly computed fields. Not computed automatically over
        % subsetting operations.
        
        % Magnitude of the velocity noise.
        n
        % n in the effective region.
        n_e
        
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
        
        % Range for computation and visualization. This region is denoted
        % the Effective Region.
        range
        % Length of range in 3 dimensions. Row vector.
        span
        
        % struct for properties of solvers.
        solver
        % struct for graphers.
        plotter
        % struct for invariant attributes of fluid.
        fluid
        % struct for proportions of unit scaling.
        scale
        % Customized struct for user to store related properties not ordained in object definition.
        data
        
        % Dimensionality of the field. 3 for 3D, 2, for 2D, minimum of 1D.
        ax
        
    end
    
    % Some rarer quantities used for convenience.
    properties (Access = private)
        % The boundary indices of the effective region given in ascending
        % order. Identical to [ones(1,3) vf.getSpan()'] except where the
        % index proceeds in the negative direction of position.
        ascLim
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
        function v = getVector(V, range)
            % Extract the vectors in a rectangular region of the grid,
            % where range = [i_0 j_0 k_0; i_f j_f k_f].
            v = squeeze(V(range(1,2): range(end,2), ...
                range(1,1): range(end,1), range(1,3): range(end,3), :));
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
            
            if ~isequal(ndims(X), 4) || ~isequal(size(X, 4), 3)
                error('Invalid grids: 4D matrices required!')
            end
            
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
            % Assume data is uniformly spaced in
            % position.
            if dims(2) > 1
                vf.xresol = X(1,2,1,1) - X(1,1,1,1);
            end
            vf.ybounds = [X(1,1,1,2) X(end,1,1,2)];
            if dims(1) > 1
                vf.yresol = X(2,1,1,2) - X(1,1,1,2);
            end
            vf.zbounds = [X(1,1,1,3) X(1,1,end,3)];
            if dims(3) > 1
                vf.zresol = X(1,1,2,3) - X(1,1,1,3);
            end
            vf.resol = [vf.xresol vf.yresol vf.zresol];
            
            % Dimensionality of our field.
            ax = sum(vf.dims > 1);
            
            vf.range = [ones(3, 1) vf.dims'];
            vf.span = (vf.range*[-1; 1])' + 1;
            vf.X_e = X;
            vf.U_e = U;
            vf.N_e = vf.N;
            
            % Ascending positions.
            vf.ascLim = [ones(3, 1) vf.getSpan()'];
            for i = 1: 3
                sp = vf.resol(i);
                if sp < 0
                    vf.ascLim(i, :) = flip(vf.ascLim(i, :));
                end
            end
            
            % Compute vorticity.
            if ax == 3
                vf.vort = vf.vorticity(0);
                vf.vort_e = vf.vort;
            end
            
            vf.initPropertyStructs()
            set(0, 'defaultTextInterpreter', 'latex');
            set(0, 'DefaultLegendInterpreter', 'latex')
        end
        
        function vfd = downsample(vf, winsize, overlap, with_noise, newXscale)
            % Perform box averaging on the effective region of the velocity
            % field and construct a downsampled field.
            %
            % newXscale is presumed to be a scalar, so that the positions
            % of the grid are in the same unit.
            
            % By default retain the coordinates.
            if ~exist('newXscale', 'var')
                newXscale = 1;
            end
            
            [Xd, Ud] = PIV_window_sim(vf.X_e, vf.U_e + with_noise*vf.N_e, ...
                winsize, overlap, newXscale);
            vfd = VelocityField(Xd, Ud);
            % Preserve certain properties.
            vfd.fluid = vf.fluid;
            vfd.scale.len = vf.scale.len / newXscale;
            vfd.derivePropertyStructs()
        end
        
        function initPropertyStructs(vf)
            % Initialize attributes encoded in struct objects for plotting,
            % computing derived quantities, and innate properties of the
            % fluid. These attributes are expected to be modified directly
            % by the user, without getter or setter functions.
            
            % Default innate attributes.
            vf.fluid.density = 1000;
            
            % Scaling factors.
            vf.scale.len = 0.001;
            
            % Properties of differentiation.
            vf.solver.diff.order = 1;
            vf.solver.diff.mode = 'central';
            
            vf.solver.ke.mode = 'direct';
            
            % Default plotter attributes.
            vf.plotter.quiverScale = 1;
            
            vf.derivePropertyStructs()
        end
        
        function derivePropertyStructs(vf)
            % Derive certain attributes from given elementary ones.
            
            vf.solver.dv = vf.scale.len^3*vf.xresol*vf.yresol*vf.zresol;
        end
            
        function setPropertyStructs(vf, fluid, solver, plotter)
            % Modify the requisite property structs.
            
            vf.fluid = fluid;
            vf.solver = solver;
            vf.plotter = plotter;
        end
        
        % Setters needed for updating dependent quantities, e.g. dv.
        
        function setRange(vf, range)
            % Set range = [i_min i_max; j_min j_max; k_min k_max] on which
            % computation and plotting are performed.
            
            % Short hand for resetting to global.
            if range == 0
                vf.setRange([1 vf.dims(2); 1 vf.dims(1); 1 vf.dims(3)])
                return
            end
            
            vf.range(1, :) = range(2, :);
            vf.range(2, :) = range(1, :);
            vf.range(3, :) = range(3, :);
            vf.span = (vf.range*[-1; 1])' + 1;
            % Ascending positions.
            vf.ascLim = [ones(3, 1) vf.getSpan()'];
            for i = 1: 3
                sp = vf.resol(i);
                if sp < 0
                    vf.ascLim(i, :) = flip(vf.ascLim(i, :));
                end
            end
            
            % Subset effective region.
            vf.X_e = vf.subsetVector(vf.X);
            vf.U_e = vf.subsetVector(vf.U);
            vf.N_e = vf.subsetVector(vf.N);
            vf.vort_e = vf.subsetVector(vf.vort);
        end
        
        function setRangePosition(vf, Xrange)
            % Variant of vf.setRange that first derives range in indices
            % from given positions.
            range(1, :) = vf.getIndex_x(Xrange(1, :));
            range(2, :) = vf.getIndex_y(Xrange(2, :));
            range(3, :) = vf.getIndex_z(Xrange(3, :));
            
            if (prod(range(:,2) - range(:,1)) < 0)
                error('Index of Subgrid proceeding in negative indices')
            end
            
            vf.setRange(range)
        end
        
        function range = getRange(vf)
            % Obtain the current effective region with indices organized in
            % x y z order.
            
            range = [vf.range(2,:); vf.range(1,:); vf.range(3,:)];
        end
        
        function dims = getDims(vf)
            % Size in x y z order of the entire grid.
            
            dims = [vf.dims(2) vf.dims(1) vf.dims(3)];
        end
        
        function dims = getSpan(vf)
            % Size in x y z order of effective region.
            
            dims = [vf.span(2) vf.span(1) vf.span(3)];
        end
        
        function v = subsetVector(vf, V)
            % Identical to VF.getVector except that 'range' is now in the
            % internal format of vf.range = [j_0 j_f; i_0 i_f; k_0 k_f].
            v = V(vf.range(1,1): vf.range(1,2), ...
                vf.range(2,1): vf.range(2,2), vf.range(3,1): vf.range(3,2), :);
        end

        %%%%%%%%%%%%%%%%%%%%% Coordinate Helpers %%%%%%%%%%%%%%%%%%%%%
        % Convert spatial coordinates to indices. The position is rounded
        % to the nearest coordinate corresponding to an index.
        
        function i = getIndex_x(vf, x)
            i = round((x - vf.xbounds(1))/vf.xresol) + 1;
            % Wrap around boundary.
            i = min(i, vf.dims(2));
            i = max(i, 1);
        end
        
        function j = getIndex_y(vf, y)
            j = round((y - vf.ybounds(1))/vf.yresol) + 1;
            % Wrap around boundaries.
            j = min(j, vf.dims(1));
            j = max(j, 1);
        end
        
        function k = getIndex_z(vf, z)
            k = round((z - vf.zbounds(1))/vf.zresol) + 1;
            % Wrap around boundary.
            k = min(k, vf.dims(3));
            k = max(k, 1);
        end
        
        function ind = getIndices(vf, pos)
            ind = [vf.getIndex_x(pos(1)), vf.getIndex_y(pos(2)), vf.getIndex_z(pos(3))];
        end
        
        function int = inBounds(vf, pos)
            int = pos(:,1)>=vf.xbounds(1) && pos(:,1)<=vf.xbounds(2) && ...
                pos(:,2)>=vf.ybounds(1) && pos(:,2)<=vf.ybounds(2) && ...
                pos(:,3)>=vf.zbounds(1) && pos(:,3)<=vf.zbounds(2);
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
            % For example, index = [0 0 k], where k is the index of the
            % plane in the appropriate dimension.
            %
            % Non-vectorized, index a 2 x 3 matrix.
            
            eq = zeros(3, 2);
            eq(:, 1) = index/sum(index);
            % Obtain a position on the plane.
            eq(:, 2) = VelocityField.getVector(vf.X, index + (index==0));
        end
        
        function eqs = getRegPlaneEqs(vf, indices)
            % Index vectors stacked as rows in a 2D array.
            
            eqs = zeros(size(indices, 1), 3, 2);
            for i = 1: size(indices, 1)
                eqs(i, :, :) = vf.getRegPlaneEq(indices(i, :));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%% Modify Velocity %%%%%%%%%%%%%%%%%%%%%%
        function addVelocity(vf, U_e)
            % In place addition of velocity field in the effective region.
            
            % Accept constant shift given as a column 3-vector.
            if isequal(size(squeeze(U_e)), [3 1])
                U_e = dimen(U_e, vf.dims);
            elseif ~isequal(size(U_e), size(vf.U_e))
                error('Mismatching Dimensions of Noise and Global Velocity')
            end
            vf.U_e = vf.U_e + U_e;
            vf.U(vf.range(1,1): vf.range(1,2), vf.range(2,1): vf.range(2,2), ...
                vf.range(3,1): vf.range(3,2), :) = vf.U_e;
        end
        
        %%%%%%%%%%%%%%%%%%%%%% Add Noise %%%%%%%%%%%%%%%%%%%%%%
        
        % Introduce species of noise to velocity. Noises generated are added to the
        % present level of noise, not replacing it.
        
        function clearNoise(vf)
            vf.N = zeros(size(vf.U));
            vf.N_e = zeros(size(vf.U_e));
        end
        
        function setNoise(vf, N_e)
            % Set noise in the current effective region.
            
            if ~isequal(size(N_e), size(vf.N_e))
                error('Mismatching Dimensions of Noise and Global Velocity')
            end
            vf.N_e = N_e;
            vf.N(vf.range(1,1): vf.range(1,2), vf.range(2,1): vf.range(2,2), ...
                vf.range(3,1): vf.range(3,2), :) = N_e;
        end
        
        function noiseMagnitude(vf)
            % Compute the magnitude of velocity in effective region and
            % globally.
            vf.n = sqrt(sum(vf.N.^2, 4));
            vf.n_e = sqrt(sum(vf.N_e.^2, 4));
        end
        
        function N_e = noise_uniform(vf, mag)
            % Add uniform noise to the effective region.
            
            % Option for adding noise with magnitude defined per position.
            if isequal(size(mag), vf.span)
                mag(:,:,:,2) = mag;
                mag(:,:,:,3) = mag(:,:,:,1);
            end
            dims = (vf.range(:,2) - vf.range(:,1) + 1)';
            N_e = (rand([dims 3])*2 - 1) .* mag/sqrt(3);
            vf.N_e = vf.N_e + N_e;
            vf.N(vf.range(1,1):vf.range(1,2), vf.range(2,1):vf.range(2,2), ...
                vf.range(3,1):vf.range(3,2), :) = vf.N_e;
        end
        
        function N_e = noise_wgn(vf, sd, snr)
            % Add white gaussian noise specified by either a standard
            % deviation or a signal-to-noise ratio to the effective region.
            
            if sd == 0
                N_e = awgn(vf.U_e, snr) - vf.U_e;
                vf.N_e = vf.N_e + N_e;
                vf.N(vf.range(1,1):vf.range(1,2), vf.range(2,1):vf.range(2,2), ...
                    vf.range(3,1):vf.range(3,2), :) = vf.N_e;
            else
                dims = (vf.range(:,2) - vf.range(:,1) + 1)';
                N_e = sd/sqrt(3)*randn([dims 3]);
                vf.N_e = vf.N_e + N_e;
                vf.N(vf.range(1,1):vf.range(1,2), vf.range(2,1):vf.range(2,2), ...
                    vf.range(3,1):vf.range(3,2), :) = vf.N_e;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%% Smoothers %%%%%%%%%%%%%%%%%%%%%%%%%
        
        function U_pre = smoothVelocity(vf, smoother)
            % Inplace smoothing of velocity field. Supported options of
            % 'smoother' are "box" and "gaussian".
            
            % Return previous velocity for convenience.
            U_pre = vf.U_e;
            
            % Cannot smooth over a unit volume.
            if max(size(vf.U_e, 1:3)) == 1
                return
            end
            
            vf.U_e(:,:,:,1) = smooth3(vf.U_e(:,:,:,1), smoother);
            vf.U_e(:,:,:,2) = smooth3(vf.U_e(:,:,:,2), smoother);
            vf.U_e(:,:,:,3) = smooth3(vf.U_e(:,:,:,3), smoother);
            % Update global velocity.
            vf.U(vf.range(1,1):vf.range(1,2), vf.range(2,1):vf.range(2,2), ...
                vf.range(3,1):vf.range(3,2), :) = vf.U_e;
        end
            
        
        function N_pre = smoothNoise(vf, smoother)
            % Alter only the noise parameter to simulate the behavior of
            % smoothing without modifying the velocity.
            
            % Return the original noise for convenience.
            N_pre = vf.N_e;
            
            % Cannot smooth over a unit volume.
            if max(size(vf.U_e, 1:3)) == 1
                return
            end
            
            vf.N_e(:,:,:,1) = smooth3(vf.U_e(:,:,:,1) + vf.N_e(:,:,:,1), smoother) - vf.U_e(:,:,:,1);
            vf.N_e(:,:,:,2) = smooth3(vf.U_e(:,:,:,2) + vf.N_e(:,:,:,2), smoother) - vf.U_e(:,:,:,2);
            vf.N_e(:,:,:,3) = smooth3(vf.U_e(:,:,:,3) + vf.N_e(:,:,:,3), smoother) - vf.U_e(:,:,:,3);
            vf.setNoise(vf.N_e)
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%% Plotters %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function plt = plotVector(vf, V, noise, title_str)
            % Make a quiver plot of the given vector field over the range
            % of interest on the current grid. 'V' can be either substted
            % already from the entire grid or be equal to 'vf.X' in
            % dimension, which will be subsetted automatically based on
            % 'vf.range'.
            
            if ~isequal(size(V, 1:3), vf.span)
                V = vf.subsetVector(V);
            end
            % Global noise allowed for automatic subsetting.
            if isequal(size(noise, 1:3), vf.dims)
                noise = vf.subsetVector(noise);
            end
            
            plt = plotVF(vf.X_e, V + noise, vf.plotter.quiverScale);
            title(title_str)
        end
        
        function plt = plotScalar(vf, S, noise, title_str)
            % Show a color map of the scalar field over the range
            % specified, defaulted to global.
            
            if ~isequal(size(S, 1:3), vf.span)
                S = vf.subsetVector(S);
            end
            % Global noise allowed for automatic subsetting.
            if isequal(size(noise, 1:3), vf.dims)
                noise = vf.subsetVector(noise);
            end
            
            % Collapse 4D matrix into a set of points.
            X = reshape(vf.X_e, [], 3);
            S = S(:);
            % Scaled size of each dot displayed. TODO: better fitting
            % expression for dot size.
            dif = sort(vf.range(:,2) - vf.range(:,1), 2);
            dot_size = 40^(sum(dif ~= 0)-1)/(dif(1)*dif(2));
            
            plt = figure;
            scatter3(X(:,1), X(:,2), X(:,3), dot_size, S + noise, 'filled');
            
            colorbar
            xlabel('$x$')
            ylabel('$y$')
            zlabel('$z$')
            title(title_str)
        end
        
        function plt = slicePlanes(vf, S, planes, noise, title_str)
            % SLICEPLANES renders continuous color plots on the planes
            % whose equations are specified. planes are in the format of
            % [[orth base]...] for all the planes to be plotted.
            
            plt = figure;
            % Subset grid.
            if ~isequal(size(S, 1:3), vf.span)
                S = vf.subsetVector(S);
            end
            % Global noise allowed for automatic subsetting.
            if isequal(size(noise, 1:3), vf.dims)
                noise = vf.subsetVector(noise);
            end
            
            for i = 1: size(planes, 1)
                eq = squeeze(planes(i, :, :));
                % Compute z value based on equation of plane.
                z = (dot(eq(:,2), eq(:,1)) - eq(1,1)*vf.X_e(:,:,1,1) - ...
                    eq(2,1)*vf.X_e(:,:,1,2)) / eq(3,1);
                % Guard against planes with unrestricted z.
                if isequal(isinf(z)+isnan(z), ones(size(z)))
                    % Regular planes.
                    if max(size(find(eq(:, 1)))) == 1
                        xslice = [];
                        yslice = [];
                        zslice = [];
                        switch find(eq(:, 1))
                            case 1
                                xslice = [eq(1, 2)];
                            case 2
                                yslice = [eq(2, 2)];
                        end
                        slice(squeeze(vf.X_e(:,:,:,1)), squeeze(vf.X_e(:,:,:,2)), ...
                            squeeze(vf.X_e(:,:,:,3)), S + noise, ...
                            xslice, yslice, zslice);
                    else
                        % Else the x (or y) value cannot be unconstrained.
                        % Derive x from y and z.
                        x = (dot(eq(:,2), eq(:,1)) - eq(2,1)*vf.X_e(:,1,:,2) - ...
                            eq(3,1)*vf.X_e(:,1,:,3)) / eq(1,1);
                        size(x)
                        size(vf.X_e(:,1,:,2))
                        size(squeeze(vf.X_e(:,1,:,3)))
                        slice(squeeze(vf.X_e(:,:,:,1)), squeeze(vf.X_e(:,:,:,2)), ...
                            squeeze(vf.X_e(:,:,:,3)), S + noise, ...
                            x, squeeze(vf.X_e(:,1,:,2)), squeeze(vf.X_e(:,1,:,3)));
                    end
     
                else
                    slice(squeeze(vf.X_e(:,:,:,1)), squeeze(vf.X_e(:,:,:,2)), ...
                        squeeze(vf.X_e(:,:,:,3)), S + noise, ...
                        squeeze(vf.X_e(:,:,1,1)), squeeze(vf.X_e(:,:,1,2)), z);
                end
                hold on
            end
            
            colorbar
            xlabel('$x$')
            ylabel('$y$')
            zlabel('$z$')
            title(title_str)
        end
        
        function plt = plotPlaneVector(vf, V, eq, noise, title_str)
            % Plots an arbitrary plane in 3D space given either three
            % non-colinear points or a normal vector + base position paris
            % formula, stored as column vectors in a matrix.
            
            if ~isequal(size(V, 1:3), vf.span)
                V = vf.subsetVector(V);
            end
            % Global noise allowed for automatic subsetting.
            if isequal(size(noise, 1:3), vf.dims)
                noise = vf.subsetVector(noise);
            end
            % Obtain a matching 4D boolean matrix indicating membership of
            % points on the plane.
            onPlane = skewPlaneMatrix(vf.X_e, eq(:,1), eq(:,2), 3);
            
            plt = plotVF(vf.X_e, (V + noise) .* onPlane, vf.plotter.quiverScale);
            title(title_str)
        end
        
        function plt = plotScalarPlaneSkewed(vf, Mag, eq, noise, title_str)
            % Derive plane equation if a normal plane is given.
            if isequal(eq(:, 2), zeros(3, 1))
                
            end
            % Boolean predicate for determining if points are on plane.
            onPlane = skewPlaneMatrix(vf.X, eq(:, 1), eq(:, 2), 1);
        end
        
        function plt = plotPlaneScalar(vf, S, range, noise, title_str)
            % Plots a scalar field S over X. This method requires that S is
            % given over the entire grid. To be fixed?
            
            % S is a 3D matrix of corresponding to positions in X.
            % Only a regular plane is currently allowed.
            % 'range' here is in the standard xyz format, of three rows
            % each giving the beginning and ending indices in that
            % dimension, not the flipped meshgrid format.
            
            % Global noise allowed for automatic subsetting.
            if isequal(size(noise, 1:3), vf.dims)
                noise = vf.subsetVector(noise);
            end
            
            % Flip x-y due to meshgrid convention.
            yrange = range(2, :);
            range(2, :) = range(1, :);
            range(1, :) = yrange;
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
            colorbar
            title(title_str)
        end
        
        function plt = isosurfaces(vf, S, isovals, noise, title_str)
            % Given an array of values of a scalar field S, their
            % isosurfaces are plotted.
            
            % Subset region of interest.
            if ~isequal(size(S, 1:3), vf.span)
                S = vf.subsetVector(S);
            end
            % Global noise allowed for automatic subsetting.
            if isequal(size(noise, 1:3), vf.dims)
                noise = vf.subsetVector(noise);
            end
            
            plt = figure;
            for isoval = isovals
                isosurface(vf.X_e(:,:,:,1), vf.X_e(:,:,:,2), vf.X_e(:,:,:,3), ...
                    S + noise, isoval)
                hold on
            end
            colorbar
            xlabel('$x$')
            ylabel('$y$')
            zlabel('$z$')
            title(title_str)
        end
        
        %%%%%%%%%%%%%%%%%%% Solvers of Derived Quantities %%%%%%%%%%%%%%%%%%
        % Some of these quantities, e.g. vorticity, is pre-computed during
        % initialization from the position and velocity field given.
        
        function vort = vorticity(vf, with_noise)
            % Computes and stores the vorticity with unit.
            % 
            % Currently just a wrapper for the built-in curl function.
            
            [vort(:,:,:,1), vort(:,:,:,2), vort(:,:,:,3)] = ...
                curl(vf.X_e(:,:,:,1), vf.X_e(:,:,:,2), vf.X_e(:,:,:,3), ...
                vf.U_e(:,:,:,1) + with_noise*vf.N_e(:,:,:,1), ...
                vf.U_e(:,:,:,2) + with_noise*vf.N_e(:,:,:,1), ...
                vf.U_e(:,:,:,3) + with_noise*vf.N_e(:,:,:,1));
        end
        
        function k = kineticEnergy(vf, with_noise)
            
            % Selected mode of computation.
            switch vf.solver.ke.mode
                case 'direct'
                    k = 1/2*vf.fluid.density * abs(vf.solver.dv) * vf.scale.len^2 * ...
                            sum((vf.U_e + with_noise*vf.N_e).^2, 'all');
            end
        end
        
        function u_mean = meanSpeed(vf, with_unit, with_noise)
            speed_e = sqrt(sum((vf.U_e + with_noise*vf.N_e).^2, 4));
            if with_unit
                u_mean = mean(speed_e, 'all') * vf.scale.len;
            else
                u_mean = mean(speed_e, 'all');
            end
        end
        
        function I = impulse(vf, with_noise, origin)
            % Computes the momentum of the fluid in the current region of
            % interest.
            %
            % 'origin' specifies the arbitrary reference point used in
            % computing impulse. If not specified, it is assumed to be the
            % origin implied in the current coordinates.
            %
            % 'I' is a column vector.
            
            if ~exist('origin', 'var')
                origin = [0 0 0]';
            end
            
            if with_noise
                I = squeeze(vf.fluid.density/2 * ...
                    sum(cross(vf.X_e - dimen(origin, vf.span), vf.vorticity(1), 4), [1 2 3]) * ...
                    vf.solver.dv*vf.scale.len);
            else
                I = squeeze(vf.fluid.density/2 * ...
                    sum(cross(vf.X_e - dimen(origin, vf.span), vf.vort_e, 4), [1 2 3]) ...
                    * vf.solver.dv*vf.scale.len);
            end
        end
        
        
        
        %%%%%%%%%%%%%%%%%%% Differential Methods %%%%%%%%%%%%%%%%%%%%%%
        
        function nder = diff1(vf, V, ind_inc, mode)
            % nder = vf.DIFF1(V, ind_inc, mode) computes the first order
            % (simple difference) numerical direction derivative of 'V' as
            % specified by the vector increment in index 'ind_inc'. The
            % parameter 'mode' specified the finite difference scheme used:
            % 'right', 'left', or 'central'. The derivative is computed in
            % the region of the volume for which the difference exists.

            % Subset region of interest.
            if ~isequal(size(V, 1:3), vf.span)
                V = vf.subsetVector(V);
            end
            % Flip meshgrid x-y.
            ind_inc2 = [ind_inc(2) ind_inc(1) ind_inc(3)];
            
            % Range of differentiable positions, from one of the bounds on the grid.
            rem = size(V, 1:3) - abs(ind_inc2);

            switch mode
                case 'right'
                    right_shifted = NaN(size(V));
                    % Range of differentiable indices.
                    idx = zeros(2, 3, 2);
                    for i = 1:3
                        if ind_inc2(i) < 0
                            idx(1, i, :) = [1-ind_inc2(i) size(V, i)];
                            idx(2, i, :) = [1 rem(i)];
                        else
                            idx(1, i, :) = [1 rem(i)];
                            idx(2, i, :) = [ind_inc2(i)+1 size(V, i)];
                        end
                    end
                    right_shifted(idx(1,1,1):idx(1,1,2), idx(1,2,1):idx(1,2,2), ...
                        idx(1,3,1):idx(1,3,2), :) = ...
                        V(idx(2,1,1):idx(2,1,2), idx(2,2,1):idx(2,2,2), ...
                        idx(2,3,1):idx(2,3,2), :);
                    dif = right_shifted - V;
                case 'left'
                    nder = -vf.diff1(V, -ind_inc, 'right');
                    return
                    
                case 'central'
                    nder = (vf.diff1(V, ind_inc, 'right') + vf.diff1(V, ind_inc, 'left')) / 2;
                    return
            end
            % Compute displacement in direction and differentiate.
            nder = dif / norm(vf.resol .* ind_inc);
        end
        
        function grad = gradient(vf, F)
            % Wrapper for 1st order central difference from Matlab.
            
            if ~isequal(size(F, 1:3), vf.span)
                error('Field matching the effective position grid expected!')
            end
            
            grad = NaN([size(F) 3]);
            % If a scalar field is given.
            if ndims(F) == 3
                [grad(:,:,:,1), grad(:,:,:,2), grad(:,:,:,3)] = ...
                    gradient(F, vf.xresol, vf.yresol, vf.zresol);
            % Otherwise a vector field.
            else
                % Last dimension is difference among components of one
                % vector, discard.
                [grad(:,:,:,:,1), grad(:,:,:,:,2), grad(:,:,:,:,3), ~] = ...
                    gradient(F, vf.xresol, vf.yresol, vf.zresol, 1);
            end
        end
        
        function jacob = jacobian(vf, V, mode)
            % Compute the jacobian matrix (or gradient, or derivative
            % matrix) of 'V' with numerical scheme given in 'mode'.
            
            % Subset region of interest.
            if ~isequal(size(V, 1:3), vf.span)
               V = vf.subsetVector(V);
            end
            
            jacob = NaN([size(V, 1:3) 3 size(V, 4)]);
            
            switch mode
                case 'unit'
                    jacob(:, :, :, 1, :) = (vf.xresol>0)*vf.diff1(V, [1 0 0], vf.solver.diff.mode);
                    jacob(:, :, :, 2, :) = (vf.yresol>0)*vf.diff1(V, [0 1 0], vf.solver.diff.mode);
                    jacob(:, :, :, 3, :) = (vf.zresol>0)*vf.diff1(V, [0 0 1], vf.solver.diff.mode);
            end
        end
        
        function div = div(vf, V)
            % Compute the divergence of the vector in the set region of
            % interest.
            %
            % Discrepant with expected results for now.
            
            % Subset region of interest.
            if ~isequal(size(V, 1:3), vf.span)
                V = vf.subsetVector(V);
            end
            div = NaN(size(V, 1:3));
            
            switch vf.solver.diff.order
                case 1
                    % Ensure diff direction in increasing x, y, z by boolean multiplication.
                    div = squeeze(vf.diff1(V(:,:,:,1), [1 0 0], vf.solver.diff.mode)*(vf.xresol>0) + ...
                        vf.diff1(V(:,:,:,2), [0 1 0], vf.solver.diff.mode)*(vf.yresol>0) + ...
                        vf.diff1(V(:,:,:,3), [0 0 1], vf.solver.diff.mode)*(vf.zresol>0));
            end
        end
        
        
       %%%%%%%%%%%%%%%%%%%%%%% Integral Methods %%%%%%%%%%%%%%%%%%%%%%%%
       
       % Scalar surface integral (scalar area element).
       function summed = intCubicSurf(vf, F)
           % Integrate on the faces of the current effective region, which
           % is a cube, the given field, vector or scalar, 'F', using a
           % scalar surface elements.
           % 
           % The indices of V are assumed to match that of X_e.
           
           if ~isequal(vf.span, size(F, 1:3))
               error('Mismatching Grids Dimensions!')
           end
           
           summed = abs(vf.yresol*vf.zresol)*(sum(F(:,vf.ascLim(1,1),:,:), 'all') + ...
               sum(F(:,vf.ascLim(1,2),:,:), 'all')) + ...
               abs(vf.xresol*vf.zresol)*(sum(F(vf.ascLim(2,1),:,:,:), 'all') + ...
               sum(F(vf.ascLim(2,2),:,:,:), 'all')) + ...
               abs(vf.xresol*vf.yresol)*(sum(F(:,:,vf.ascLim(3,1),:), 'all') + ...
               sum(F(:,:,vf.ascLim(3,2),:), 'all'));
           summed = squeeze(summed);
       end
       
       function vec = intCubicSurf_vec(vf, S)
           % Integrate on the faces of the current effective region, which
           % is a cube, the six vector surface elements multiplied by the
           % given scalar field.
           % 
           % The indices of V are assumed to match that of X_e.
           
           if ~isequal(vf.span, size(S, 1:3))
               error('Mismatching Grids Dimensions!')
           end
           
           vec = abs(vf.yresol*vf.zresol)*([-1 0 0]'*sum(S(:,vf.ascLim(1,1),:), 'all') + ...
               [1 0 0]'*sum(S(:,vf.ascLim(1,2),:), 'all')) + ...
               abs(vf.xresol*vf.zresol)*([0 -1 0]'*sum(S(vf.ascLim(2,1),:,:), 'all') + ...
               [0 1 0]'*sum(S(vf.ascLim(2,2),:,:), 'all')) + ...
               abs(vf.xresol*vf.yresol)*([0 0 -1]'*sum(S(:,:,vf.ascLim(3,1)), 'all') + ...
               [0 0 1]'*sum(S(:,:,vf.ascLim(3,2)), 'all'));
       end
       
       function flux = intCubicSurf_flux(vf, V)
           % Compute the flux on the six faces of the cube as defined by
           % the current effective region by the given vector field 'V'.
           % The indices of 'V' are assumed to match that of X_e.
           
           if ~isequal(vf.span, size(V, 1:3))
               error('Mismatching Grids Dimensions!')
           end
           
           % Left face.
           flux_left = sum(-V(:,vf.ascLim(1,1),:,1)*abs(vf.yresol)*abs(vf.zresol), 'all');
           % Right face.
           flux_right = sum(V(:,vf.ascLim(1,2),:,1)*abs(vf.yresol)*abs(vf.zresol), 'all');
           % Bottom face.
           flux_bottom = sum(-V(vf.ascLim(2,1),:,:,2)*abs(vf.xresol)*abs(vf.zresol), 'all');
           % Top face.
           flux_top = sum(V(vf.ascLim(2,2),:,:,2)*abs(vf.xresol)*abs(vf.zresol), 'all');
           % Back face.
           flux_back = sum(-V(:,:,vf.ascLim(3,1),3)*abs(vf.xresol)*abs(vf.yresol), 'all');
           % Front face.
           flux_front = sum(V(:,:,vf.ascLim(3,2),3)*abs(vf.xresol)*abs(vf.yresol), 'all');
           
           flux = flux_left + flux_right + flux_top + flux_bottom + flux_front + flux_back;
       end
       
       function crs = intCubicSurf_cross(vf, V)
           % Compute the integral of the cross product V x n, where the
           % latter is the normal vector on one of the six faces of the
           % effective cube. 'V' is expected as a standard 4D array, whose
           % indices are assumed to match that of X_e.
           
           if ~isequal(vf.span, size(V, 1:3))
               error('Mismatching Grids Dimensions!')
           end
           
           % Left face.
           crs_pd = zeros(vf.span(1), vf.span(3), 3);
           crs_pd(:,:,2) = -squeeze(V(:, vf.ascLim(1,1), :, 3));
           crs_pd(:,:,3) = squeeze(V(:, vf.ascLim(1,1), :, 2));
           left = abs(vf.yresol)*abs(vf.zresol)*sum(crs_pd, [1 2]);
           % Right face.
           crs_pd = zeros(vf.span(1), vf.span(3), 3);
           crs_pd(:,:,2) = squeeze(V(:, vf.ascLim(1,2), :, 3));
           crs_pd(:,:,3) = -squeeze(V(:, vf.ascLim(1,2), :, 2));
           right = abs(vf.yresol)*abs(vf.zresol)*sum(crs_pd, [1 2]);
           % Bottom face.
           crs_pd = zeros(vf.span(2), vf.span(3), 3);
           crs_pd(:,:,1) = squeeze(V(vf.ascLim(2,1), :, :, 3));
           crs_pd(:,:,3) = -squeeze(V(vf.ascLim(2,1), :, :, 1));
           bottom = abs(vf.xresol)*abs(vf.zresol)*sum(crs_pd, [1 2]);
           % Top face.
           crs_pd = zeros(vf.span(2), vf.span(3), 3);
           crs_pd(:,:,1) = -squeeze(V(vf.ascLim(2,2), :, :, 3));
           crs_pd(:,:,3) = squeeze(V(vf.ascLim(2,2), :, :, 1));
           top = abs(vf.xresol)*abs(vf.zresol)*sum(crs_pd, [1 2]);
           % Back face.
           crs_pd = zeros(vf.span(1), vf.span(2), 3);
           crs_pd(:,:,1) = -squeeze(V(:, :, vf.ascLim(3,1), 2));
           crs_pd(:,:,2) = squeeze(V(:, :, vf.ascLim(3,1), 1));
           back = abs(vf.xresol)*abs(vf.yresol)*sum(crs_pd, [1 2]);
           % Front face.
           crs_pd = zeros(vf.span(1), vf.span(2), 3);
           crs_pd(:,:,1) = squeeze(V(:, :, vf.ascLim(3,2), 2));
           crs_pd(:,:,2) = -squeeze(V(:, :, vf.ascLim(3,2), 1));
           front = abs(vf.xresol)*abs(vf.yresol)*sum(crs_pd, [1 2]);
           
           crs = squeeze(left + right + bottom + top + back + front);
       end
        
    end
end