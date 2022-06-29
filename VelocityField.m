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
        vort_e = NaN
        
        % Grid dimension. Note that the indices on grid, due to meshgrid
        % convension, is given as (y, x, z), where x represent the number
        % of distinct values of x, etc.. Row vector.
        dims
        % Optional 4D matrix of noise in velocity to be superposed with U.
        N
        % Noise in the region of interest.
        N_e
        
        % Commonly computed fields. Not computed automatically over
        % subsetting operations.
        
        % Magnitude of the velocity noise.
        n
        % n in the effective region.
        n_e
        
        % Derived fields.
        
        % Lower and upper bounds of x values.
        xbounds
        % Uniform spacing in x.
        xsp
        % Lower and upper bounds of y values.
        ybounds
        % Uniform spacing y.
        ysp
        % Lower and upper bounds of z values.
        zbounds
        % Uniform spacing in z.
        zsp
        % Array of spacings.
        sps
        % Collection of bounds, parallels vf.dims, though natural in x-y
        % ordering, which is flipped for vf.dims.
        bounds
        
        % Range for computation and visualization. This region is denoted
        % the Effective Region.
        range
        % Length of range in 3 dimensions. Row vector.
        span
        
        % Corresponding limits of position of the current effective region,
        % stored however in natural x-y-z order.
        lims
        
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
        
        % Font size for plotting.
        fontsize = 10;
    end
    
    % Some rarer quantities used for convenience.
    properties (Access = private)
        % The boundary indices of the effective region given in ascending
        % order. Identical to [ones(1,3) vf.getSpan()'] except where the
        % index proceeds in the negative direction of position.
        ascLim
    end
    
    properties (Constant)
        % For printing given dimension as number.
        dim_str = ['x', 'y', 'z'];
        % For converting between natural x, y, z dimensional ordering to
        % meshgrid y, x, z ordering.
        dim_flip = [2 1 3];
        % Flag for loading differentiation module.
        load_diff = VelocityField.load_finite_diffs();
        % Flag for setting latex font for plotting.
        set_latex = VelocityField.set_latex_fonts();
    end
    
    methods(Static)
        
        %%%%%%%%%%%%%%%%%%% Adapters of Recorded Data %%%%%%%%%%%%%%%%%%%
        function vf = importCmps(x, y, z, u, v, w, minimal)
            % Import velocity field by components of position and velocity,
            % as 3D arrays.
           
            if ~exist('minimal', 'var')
                minimal = 0;
            end
            
            dims = [size(x) 3];  
            % Pack positions and velocities compactly as 3-vectors in extra dimension.
            X = zeros(dims);
            X(:,:,:,1) = x;
            X(:,:,:,2) = y;
            X(:,:,:,3) = z;
            U = zeros(dims);
            U(:,:,:,1) = u;
            U(:,:,:,2) = v;
            U(:,:,:,3) = w;
            
            vf = VelocityField(X, U, minimal);
        end
        
        %%%%%%% Loading pre-computed differentiation formulas %%%%%%%
        function l = load_finite_diffs()
            % Load pre-computed finite difference formulas as a global
            % variable 'PDF'. This should not be replaced for centricMatrix
            % performs look-up on it.
            
            global PDF
            global rootFolder
            load(strcat(rootFolder, '\diff\finite-differences.mat'), 'PDF')
            l = 1;
        end
        
        function l = set_latex_fonts()
            % Set default font for plotting as latex.
            set(0, 'defaultTextInterpreter', 'latex');
            set(0, 'DefaultLegendInterpreter', 'latex')
            l = 1;
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
        
        %%%%%%%%%%%%%%%%%%% 4D Grid Arithmetic %%%%%%%%%%%%%%%%%%%%%%
        function V = operate3Vector(V, v, op)
            % Perform binary operation defined by 'op' on each 3-vector of
            % the 4D matrix 'V' with the given 3-vector 'v'.
            
            if ndims(V) ~= 4 || size(V, 4) ~= 3
                error('Proper 4D matrix expected!')
            end
            V(:,:,:,1) = op(V(:,:,:,1), v(1));
            V(:,:,:,2) = op(V(:,:,:,2), v(2));
            V(:,:,:,3) = op(V(:,:,:,3), v(3));
        end
        
        function V = add3Vector(V, v)
            % Add the vector to 'v' to every vector on the 3D grid of the
            % 4D array 'V'.
            
            V = VelocityField.operate3Vector(V, v, @plus);
        end
        
        function V = subtract3Vector(V, v)
            % Subtract the vector 'v' from every vector on the 3D grid of the
            % 4D array 'V', whose 4th dimension contains the vectors.
            
            V = VelocityField.operate3Vector(V, v, @minus);
        end
        
        function D = dot3Vector(V, v)
            % Take the dot product of every vector on the grid with the
            % given vector 'v'. The result is a 3D array.
            
            D = VelocityField.operate3Vector(V, v, @times);
            D = sum(D, 4);
        end
        
        function V = operateGridwise(V, r, op)
            % Perform position-wise binary operation defined by 'op' on each component
            % of the 3-vector of the 4D matrix 'V' and the corresponding scalar from
            % 'R'.
            
            V = zeros(size(V));
            if ndims(V) ~= 4 || size(V, 4) ~= 3
                error('Proper 4D matrix expected!')
            end
            V(:,:,:,1) = op(V(:,:,:,1), r);
            V(:,:,:,2) = op(V(:,:,:,2), r);
            V(:,:,:,3) = op(V(:,:,:,3), r);
        end
        
        %%%%%%%%%%%%%%%%%% Plotting Helpers %%%%%%%%%%%%%%%%%%%
        function on = skewPlaneMatrix(X, orth, base, sz)
            % SKEWPLANEMATRIX(X, orth, base) returns a matrix corresponding to the
            % positions in X, indicating whether positions therein are on the plane
            % defined by the orthogonal vector and the base position given.
            %
            % sz is the size of the last dimension used for element-wise multiplication.
            % For velocity, sz = 4; for a scalar field, say error, sz = 1.
            
            % Tolerance for rounding error.
            tol = 1e-4;
            
            % Upper-dimensionalize base and normal vectors for subtraction from 4D X.
            dims = size(X);
            dims = dims(1:3);
            
            onplane = abs(VelocityField.dot3Vector(VelocityField.subtract3Vector(X, base), ...
                orth)) < tol;
            % Expand this logical 3D matrix indicating membership on plane to
            % the size of the last dimension for direct point-wise multiplication.
            on = zeros([dims sz]);
            for i = 1: sz
                on(:,:,:,i) = onplane;
            end
        end
        
        function [orth, base, onPlane] = getPlaneEq(x)
            % [orth, base, onPlane] = getPlaneEq(X)
            %
            % Given three positions stored as columns in the x matrix, compute a
            % representation for a plane: the unit normal vector, the base position,
            % and a boolean handle that decides whether a given point is on the plane.
            
            % Pick first position as base and use displacements from it to compute a
            % unit normal vector.
            orth = cross(x(:,2) - x(:,1), x(:,3) - x(:,1));
            orth = orth / norm(orth);
            
            base = x(:,1);
            
            % Handle applicable exclusively to 4D matrix.
            onPlane = @(X) VelocityField.skewPlaneMatrix(X, orth, base);
        end
        
        %%%%%%%%%%%%%%%%%% Time-resolved Quantities %%%%%%%%%%%%%%%%%%%
        
        function ds = scalarTimeDeriv(s, dt, diff_order, err_order, noise)
            
            if ~isvector(s)
                error('A scalar array expected!')
            end
            % Ensure shape of vector.
            s = reshape(s, [], 1);
            
            if ~exist('noise', 'var')
                noise = 0;
            % Let pass for a uniform noise.
            elseif isequal(size(noise), [1 1])
            elseif ~isvector(noise) || ~isequal(length(noise), length(s))
                error('Noise not matching original vector in dimension!')
            else
                noise = reshape(noise, [], 1);
            end
            
            ds = dt^(-diff_order)*centricMatrix(size(s, 1), ...
                diff_order, err_order) * (s+noise);
        end
        
        function dv = vectorTimeDeriv(v, dt, diff_order, err_order, noise)
            % Assume the vectors are stored in the column of the matrix
            % 'v'.
            
            if ~exist('noise', 'var')
                noise = zeros(size(v));
            % Let pass for uniform noises
            elseif isequal(size(noise), [1 1])
                noise = repmat(noise, size(v, 1), size(v, 2));
            elseif isequal(size(noise), [3 1])
                noise = repmat(noise, 1, size(v, 2));
            elseif ~isequal(size(noise), size(v))
                error('Noise not matching original matrix in dimension!')
            end
            
            dv = zeros(size(v));
            for d = 1: size(v, 1)
                dv(d, :) = VelocityField.scalarTimeDeriv(v(d, :), dt, diff_order, ...
                    err_order, noise(d, :));
            end
        end
        
        function [dI_dt, I] = impulse_time_deriv(vfs, origins, dt, ...
                with_noise, smoother, span)
           % Assume regions are properly subsetted in the given array of
           % velocity fields traced over time.
           
           I = zeros(3, length(vfs));
           for i = 1: length(vfs)
               I(:,i) = vfs{i}.impulse(origins(:,i), with_noise);
           end
           
           % Smooth impulse dimensionally if demanded.
           if exist('smoother', 'var')
               for d = 1: 3
                   I(d,:) = smooth(I(d,:), span, smoother);
               end
           end
           
           % Compute derivative.
           dI_dt = zeros(3, length(vfs));
           D = centricMatrix(length(vfs), 1, vfs{1}.solver.diff.err_order) / dt;
           for d = 1: 3
               dI_dt(d,:) = D * I(d,:)';
           end
        end
        
    end
    
    methods
        
        function vf = VelocityField(X, U, minimal)
            % vf = VelocityField(X, U, minimal)
            % Accepts valid 4D matrices of position and velocity whose
            % first three dimensions conform to that of a meshgrid and last
            % dimension corresponds to 3-vectors stored in the [x y z]'
            % format.
            %
            % The flag 'minimal' when set true omits the computation of
            % derived quantities from the velocity field, currently only
            % the vorticity field.
            
            if ~isequal(ndims(X), 4) || ~isequal(size(X, 4), 3)
                error('Invalid grids: 4D matrices required!')
            end
            
            vf.X = X;
            vf.U = U;
            
            if any(isnan(U), 'all')
                warning(strcat("Some velocity vector on grid is NaN.", ...
                    " Treatment for computing differet quantities may vary!"))
            end
            
            if ~isequal(size(X), size(U))
                error('Mismatching grid dimensions for position and velocity data')
            end
            dims = size(X);
            vf.dims = dims(1:3);
            % Init null additional noise.
            vf.N = zeros(size(U));
            
            vf.xbounds = [X(1,1,1,1) X(1,end,1,1)];
            % Assume data is uniformly spaced in position.
            % Allowing for lower dimensional grid by defaulting sp to 0.
            vf.xsp = 0;
            if dims(2) > 1
                vf.xsp = X(1,2,1,1) - X(1,1,1,1);
            end
            vf.ybounds = [X(1,1,1,2) X(end,1,1,2)];
            vf.ysp = 0;
            if dims(1) > 1
                vf.ysp = X(2,1,1,2) - X(1,1,1,2);
            end
            vf.zbounds = [X(1,1,1,3) X(1,1,end,3)];
            vf.zsp = 0;
            if dims(3) > 1
                vf.zsp = X(1,1,2,3) - X(1,1,1,3);
            end
            vf.sps = [vf.xsp vf.ysp vf.zsp];
            vf.bounds = [vf.xbounds; vf.ybounds; vf.zbounds];
            vf.lims = vf.bounds;
            
            % Dimensionality of our field.
            vf.ax = sum(vf.dims > 1);
            
            vf.range = [ones(3, 1) vf.dims'];
            vf.span = (vf.range*[-1; 1])' + 1;
            vf.X_e = X;
            vf.U_e = U;
            vf.N_e = vf.N;
            
            % Ascending positions.
            vf.ascLim = [ones(3, 1) vf.getSpan()'];
            for i = 1: 3
                sp = vf.sps(i);
                if sp < 0
                    vf.ascLim(i, :) = flip(vf.ascLim(i, :));
                end
            end
            
            % Initialize innate quantities and customized settings.
            vf.initPropertyStructs()
            
            % Derive quantities if not proscribed.
            if exist('minimal', 'var') && minimal
                vf.vort_e = NaN;
            else
                vf.deriveQuantities();
            end
            
            vf.vort = vf.vort_e;
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
            % Preserve effective region.
            vfd.setRangePosition(vf.getRangePosition())
        end
        
        function deriveQuantities(vf)
            % Derive physical quantities of the velocity field in the
            % effective region. The global fields, when these exist, are
            % not updated.
            
            vf.vort_e = vf.vorticity(0);
%             if vf.ax == 3
%                 
%             elseif vf.ax ==2
%                 
%             end
        end
        
        function updateFields(vf)
            % Succeeding a modification of the effective velocity, update
            % all fields. This is expensive for vorticity computation.
            
            % Update velocity.
            vf.U(vf.range(1,1): vf.range(1,2), vf.range(2,1): vf.range(2,2), ...
                vf.range(3,1): vf.range(3,2), :) = vf.U_e;
            vf.vort_e = vf.vorticity(0);
        end
        
        function updateGlobal(vf)
            % Update the global fields, excepting the velocity field, based
            % on the field in the current effective region. This is,
            % updatting the field data structures with suffix 'e' to the
            % corresponding portions in their global counterparts.
            %
            % Given that this handle updates all the data structures of
            % fields supported in this class, its use is expensive and
            % dissuaded unless modification is resulted in all such fields.
            
            % Update velocity.
            vf.U(vf.range(1,1): vf.range(1,2), vf.range(2,1): vf.range(2,2), ...
                vf.range(3,1): vf.range(3,2), :) = vf.U_e;
            % Update noise.
            vf.N(vf.range(1,1): vf.range(1,2), vf.range(2,1): vf.range(2,2), ...
                vf.range(3,1): vf.range(3,2), :) = vf.N_e;
            % Update physical quantities.
            vf.vort(vf.range(1,1): vf.range(1,2), vf.range(2,1): vf.range(2,2), ...
                vf.range(3,1): vf.range(3,2), :) = vf.vort_e;
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
            vf.solver.diff.err_order = 2;
            vf.solver.diff.mode = 'centric';
            
            % Default plotter attributes.
            vf.plotter.quiverScale = 1;
            
            vf.derivePropertyStructs()
        end
        
        %%%%%%%%%%%%%%%%%%%%% Set Physical Units %%%%%%%%%%%%%%%%%%%%%%
        function setLengthScale(vf, toMeter)
            vf.scale.len = toMeter;
            vf.derivePropertyStructs()
        end
        
        function derivePropertyStructs(vf)
            % Derive certain attributes from given elementary ones.
            
            % Allow lower dimensional element, area or length, for which
            % the unit is assumed to be the same as the largest non-trivial
            % grid length.
            sps = vf.sps(vf.sps~=0);
            vf.solver.dv = vf.scale.len^3*prod(sps) * max(sps)^(3-vf.ax);
        end
            
        function setPropertyStructs(vf, fluid, solver, plotter)
            % Modify the requisite property structs.
            
            vf.fluid = fluid;
            vf.solver = solver;
            vf.plotter = plotter;
        end
        
        % Setters needed for updating dependent quantities, e.g. dv.
        
        function setRange(vf, range, minimal)
            % Set range = [i_min i_max; j_min j_max; k_min k_max] on which
            % computation and plotting are performed.
            %
            % Note that differentiated quantities, such as vorticity, are
            % re-calculated instead of simply subsetted.
            
            % Short hand for resetting to global.
            if range == 0
                vf.setRange([1 vf.dims(2); 1 vf.dims(1); 1 vf.dims(3)])
                return
            end
            
            vf.range(1, :) = range(2, :);
            vf.range(2, :) = range(1, :);
            vf.range(3, :) = range(3, :);
            vf.span = (vf.range*[-1; 1])' + 1;
            vf.lims = vf.getX(range);
            
            % Ascending positions.
            vf.ascLim = [ones(3, 1) vf.getSpan()'];
            for i = 1: 3
                sp = vf.sps(i);
                if sp < 0
                    vf.ascLim(i, :) = flip(vf.ascLim(i, :));
                end
            end
            
            % Subset effective region.
            vf.X_e = vf.subsetField(vf.X);
            vf.U_e = vf.subsetField(vf.U);
            vf.N_e = vf.subsetField(vf.N);
            
            % Recompute vorticity since the boundary values are different.
            if exist('minimal', 'var') && minimal == 1
                vf.anullQuantities()
            else
                vf.deriveQuantities()
                vf.vort_e = vf.subsetField(vf.vort);
            end
        end
        
        function setRangePosition(vf, Xrange, minimal)
            % Variant of vf.setRange that first derives range in indices
            % from given positions. 'Xrange' can be given in either
            % ascending or descending order per row.
            
            range(1, :) = vf.getIndex_x(Xrange(1, :));
            range(2, :) = vf.getIndex_y(Xrange(2, :));
            range(3, :) = vf.getIndex_z(Xrange(3, :));
            range = sort(range, 2);
            
            if prod(range(:,2) - range(:,1)) == 0
                error('Not a 3D subgrid!')
            end
            
            % If physical quantities should not be automatically derived.
            if ~exist('minimal', 'var')
                minimal = 0;
            end
            vf.setRange(range, minimal)
        end
        
        function range = getRange(vf)
            % Obtain the current effective region with indices organized in
            % x y z order.
            
            range = [vf.range(2,:); vf.range(1,:); vf.range(3,:)];
        end
        
        function Xrange = getRangePosition(vf)
            Xrange = [vf.X_e(1,1,1,1) vf.X_e(1,end,1,1); ...
                    vf.X_e(1,1,1,2) vf.X_e(end,1,1,2); ...
                    vf.X_e(1,1,1,3) vf.X_e(1,1,end,3)];
        end
        
        function dims = getDims(vf)
            % Size in x y z order of the entire grid.
            
            dims = [vf.dims(2) vf.dims(1) vf.dims(3)];
        end
        
        function dims = getSpan(vf)
            % Size in x y z order of effective region.
            
            dims = [vf.span(2) vf.span(1) vf.span(3)];
        end
        
        function v = subsetField(vf, V)
            % Identical to VF.getVector except that 'range' is now in the
            % internal format of vf.range = [j_0 j_f; i_0 i_f; k_0 k_f].
            % Both scalar and vector fields allowed.
            v = V(vf.range(1,1): vf.range(1,2), ...
                vf.range(2,1): vf.range(2,2), vf.range(3,1): vf.range(3,2), :);
        end

        %%%%%%%%%%%%%%%%%%%%% Coordinate Helpers %%%%%%%%%%%%%%%%%%%%%
        % Convert spatial coordinates to indices. The position is rounded
        % to the nearest coordinate corresponding to an index.
        
        function v = getVectorat(vf, V, pos)
            v = VelocityField.getVector(V, vf.getIndices(pos));
        end
        
        function i = getIndex_x(vf, x)
            i = round((x - vf.xbounds(1))/vf.xsp) + 1;
            % Wrap around boundary.
            i = min(i, vf.dims(2));
            i = max(i, 1);
        end
        
        function j = getIndex_y(vf, y)
            j = round((y - vf.ybounds(1))/vf.ysp) + 1;
            % Wrap around boundaries.
            j = min(j, vf.dims(1));
            j = max(j, 1);
        end
        
        function k = getIndex_z(vf, z)
            k = round((z - vf.zbounds(1))/vf.zsp) + 1;
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
            x = vf.xbounds(1) + (i-1)*vf.xsp;
        end
        
        function y = get_y(vf, j)
            y = vf.ybounds(1) + (j-1)*vf.ysp;
        end
        
        function z = get_z(vf, k)
            z = vf.zbounds(1) + (k-1)*vf.zsp;
        end
        
        function X = getX(vf, index)
            X = [vf.get_x(index(1,:)); vf.get_y(index(2,:));  vf.get_z(index(3,:))];
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
                vf.U_e = VelocityField.add3Vector(vf.U_e, U_e);
            elseif ~isequal(size(U_e), size(vf.U_e))
                error('Mismatching Dimensions of Noise and Global Velocity')
            else
                vf.U_e = vf.U_e + U_e;
            end
            % Update velocity.
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
        
        function setFontSize(vf, fsize)
            vf.fontsize = fsize;
        end
        
        function plt = plotVector(vf, V, noise, title_str, varargin)
            % Make a quiver plot of the given vector field over the range
            % of interest on the current grid. 'V' can be either substted
            % already from the entire grid or be equal to 'vf.X' in
            % dimension, which will be subsetted automatically based on
            % 'vf.range'. 'varargin' takes in typical name value pairs of
            % quiver3 options as passed into the Matlab quiver3 function.
            %
            % Note that if the plot is desired to be generated on a new
            % figure, that must be manually created before invoking this
            % function.
            
            
            if ~isequal(size(V, 1:3), vf.span)
                V = vf.subsetField(V);
            end
            % Global noise allowed for automatic subsetting.
            if isequal(size(noise, 1:3), vf.dims)
                noise = vf.subsetField(noise);
            end
            
            plt = plotVF(vf.X_e, V + noise, vf.plotter.quiverScale, vf.range, varargin{:});
            title(title_str, 'FontSize', vf.fontsize)
        end
        
        function plt = plotScalar(vf, S, noise, title_str)
            % Show a color map of the scalar field over the range
            % specified, defaulted to global.
            
            if ~isequal(size(S, 1:3), vf.span)
                S = vf.subsetField(S);
            end
            % Global noise allowed for automatic subsetting.
            if isequal(size(noise, 1:3), vf.dims)
                noise = vf.subsetField(noise);
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
            xlabel('$x$', 'FontSize', vf.fontsize)
            ylabel('$y$', 'FontSize', vf.fontsize)
            zlabel('$z$', 'FontSize', vf.fontsize)
            title(title_str, 'FontSize', vf.fontsize)
        end
        
        function plt = slicePlanes(vf, S, planes, noise, title_str)
            % SLICEPLANES renders continuous color plots on the planes
            % whose equations are specified. planes are in the format of
            % [[orth base]...] for all the planes to be plotted.
            
            plt = figure;
            % Subset grid.
            if ~isequal(size(S, 1:3), vf.span)
                S = vf.subsetField(S);
            end
            % Global noise allowed for automatic subsetting.
            if isequal(size(noise, 1:3), vf.dims)
                noise = vf.subsetField(noise);
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
            % Plots an arbitrary plane in 3D space given as a 
            % (normal vector, base position) pair, stored as column vectors
            % in the matrix eq. To obtain such a definition of a plane from
            % three non-colinear points, use vf.getPlaneEq().
            
            if ~isequal(size(V, 1:3), vf.span)
                V = vf.subsetField(V);
            end
            % Global noise allowed for automatic subsetting.
            if isequal(size(noise, 1:3), vf.dims)
                noise = vf.subsetField(noise);
            end
            % Obtain a matching 4D boolean matrix indicating membership of
            % points on the plane.
            onPlane = VelocityField.skewPlaneMatrix(vf.X_e, eq(:,1), eq(:,2), 3);
            
            plt = plotVF(vf.X_e, (V + noise) .* onPlane, vf.plotter.quiverScale);
            title(title_str)
        end
        
        function plt = plotScalarPlaneSkewed(vf, Mag, eq, noise, title_str)
            % Derive plane equation if a normal plane is given.
            if isequal(eq(:, 2), zeros(3, 1))
                
            end
            % Boolean predicate for determining if points are on plane.
            onPlane = VelocityField.skewPlaneMatrix(vf.X, eq(:, 1), eq(:, 2), 1);
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
                noise = vf.subsetField(noise);
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
                S = vf.subsetField(S);
            end
            % Global noise allowed for automatic subsetting.
            if isequal(size(noise, 1:3), vf.dims)
                noise = vf.subsetField(noise);
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
            % Compute the vorticity field with unit. NaN values in velocity
            % will result in NaN in vorticity, which are not replaced.
            
            vort = vf.curl(vf.U_e + with_noise*vf.N_e);
        end
        
        function K = kineticEnergy(vf, with_noise)
            % Compute the kinetic energy of the effective region with
            % definitional formula. NaN values of velocity are ignored in
            % sum.
            
            K = 1/2*vf.fluid.density * abs(vf.solver.dv) * vf.scale.len^2 * ...
                            sum((vf.U_e + with_noise*vf.N_e).^2, 'all', 'omitnan');
        end
        
        function u_mean = meanSpeed(vf, with_unit, with_noise)
            speed_e = sqrt(sum((vf.U_e + with_noise*vf.N_e).^2, 4));
            if with_unit
                u_mean = mean(speed_e, 'all') * vf.scale.len;
            else
                u_mean = mean(speed_e, 'all');
            end
        end
        
        function I = impulse(vf, origin, with_noise)
            % Compute the momentum of the fluid in the current region of
            % interest. NaN entires in vorticity field is ignored in total
            % sum over effective region.
            %
            % 'origin' specifies the arbitrary reference point used in
            % computing impulse. If not specified, it is assumed to be the
            % origin implied in the current coordinates.
            %
            % 'I' returned is a column vector.
            
            if ~isequal(size(origin), [3 1])
                error('Invalid origin')
            end
            
            if ~exist('origin', 'var')
                origin = [0 0 0]';
            end
            
            if with_noise
                I = squeeze(vf.fluid.density/2 * ...
                    sum(cross(VelocityField.subtract3Vector(vf.X_e, origin), ...
                        vf.vorticity(1), 4), [1 2 3], 'omitnan') * ...
                    vf.solver.dv*vf.scale.len);
            else
                I = squeeze(vf.fluid.density/2 * ...
                    sum(cross(VelocityField.subtract3Vector(vf.X_e, origin), ...
                        vf.vort_e, 4), [1 2 3], 'omitnan') * ...
                    vf.solver.dv*vf.scale.len);
            end
        end
        
        function flux = impulse_flux(vf, origin, with_noise, out)
            % computes the impulse flux on the surfaces of the volume
            
            % 'origin' specifies the arbitrary reference point used in
            % computing impulse. If not specified, it is assumed to be the
            % origin implied in the current coordinates.
            
            % 'out' specifies the format for returning quantities:
            % 'sum' returns the total impulse flux over all surfaces
            % of the volume
            % 'faces' returns the impulse flux on each face in the
            % order: [LEFT RIGHT BOTTOM TOP BACK FRONT]
            
            if ~isequal(size(origin), [3 1])
                error('Invalid origin')
            end
            
            if ~exist('origin', 'var')
                origin = [0 0 0]';
            end
            
            if ~exist('out', 'var')
                out = 'sum';
            end
            
            if with_noise
                I = vf.fluid.density/2 * ...
                    cross(VelocityField.subtract3Vector(vf.X_e, origin), ...
                    vf.vorticity(1), 4)*vf.scale.len;
           
            else
                I = vf.fluid.density/2 * ...
                    cross(VelocityField.subtract3Vector(vf.X_e, origin), ...
                    vf.vort_e, 4)*vf.scale.len;
            end
            
            fluxX = (vf.scale.len)^3*intCubicSurf_flux(vf, vf.U_e.*I(:,:,:,1), out);
            fluxY = (vf.scale.len)^3*vf.scale.len*intCubicSurf_flux(vf, vf.U_e.*I(:,:,:,2), out);
            fluxZ = (vf.scale.len)^3*vf.scale.len*intCubicSurf_flux(vf, vf.U_e.*I(:,:,:,3), out);
            
            switch out
                case 'sum'
                    flux = [sum(fluxX,'omitnan'); sum(fluxY,'omitnan'); sum(fluxZ,'omitnan')];
                case 'faces'
                    flux = [fluxX; fluxY; fluxZ];
            end
        end
        
        function flux = kineticEnergy_flux(vf, with_noise, out)
            % computes the energy flux on the surfaces of the volume
                       
            % 'out' specifies the format for returning quantities:
            % 'sum' returns the total impulse flux over all surfaces
            % of the volume
            % 'faces' returns the impulse flux on each face in the
            % order: [LEFT RIGHT BOTTOM TOP BACK FRONT]
            
            K = sum(1/2*vf.fluid.density * vf.scale.len^2 * ...
                (vf.U_e + with_noise*vf.N_e).^2, 4, 'omitnan');
            
            flux = (vf.scale.len)^3*intCubicSurf_flux(vf, vf.U_e.*K, out);

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
                V = vf.subsetField(V);
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
            nder = dif / norm(vf.sps .* ind_inc);
        end
        
        %%%%%%%%%%%%% Spatial Differentiation %%%%%%%%%%%%%
        
        function dF = diff(vf, F, dim, diff_order)
            % Differentiate the given scalar or vector field 'F' with
            % respect to change in 'dim' dimension, using polynomial
            % interpolant of order 'n'.
            
            % Order of accuracy.
            err_order = vf.solver.diff.err_order;
            % Differentiation matrix to be applied on each edge.
            D = 1/vf.sps(dim)^diff_order * ...
                centricMatrix(size(F, vf.dim_flip(dim)), diff_order, err_order);
            
            dF = zeros(size(F));
            
            % Differentiate each dimension.
            for d = 1: size(F, 4)
                switch dim
                    case 1
                        for j = 1: size(F, 1)
                            for k = 1: size(F, 3)
                                dF(j,:,k,d) = D*reshape(F(j,:,k,d), [], 1);
                            end
                        end
                    case 2
                        for i = 1: size(F, 2)
                            for k = 1: size(F, 3)
                                dF(:,i,k,d) = D*reshape(F(:,i,k,d), [], 1);
                            end
                        end
                    case 3
                        for j = 1: size(F, 1)
                            for i = 1: size(F, 2)
                                dF(j,i,:,d) = D*reshape(F(j,i,:,d), [], 1);
                            end
                        end
                end
            end
        end
        
        function grad = gradient(vf, F)
            % Spatial gradient of the given field 'F' with grid spacing
            % corresponding to the measurement grid. If a 3D field is
            % given, the vector gradient is computed; for a 4D field, the
            % Jacobian matrix is returned.
            
            if ~isequal(size(F, 1:3), vf.span)
                error('Field matching the effective position grid expected!')
            end
            
            % By default do not omit undefined derivatives.
            grad = NaN([size(F) 3]);
            
            % If a scalar field is given.
            if ndims(F) == 3
                for j = 1: 3
                    grad(:,:,:,j) = vf.diff(F, j, 1);
                end
            % Otherwise a vector field.
            else
                for i = 1: size(F, 4)
                    for j = 1: 3
                        grad(:,:,:,i,j) = vf.diff(F(:,:,:,i), j, 1);
                    end
                end
            end
        end
        
        
        function div = div(vf, V)
            % Compute the divergence of the vector in the set region of
            % interest.
            %
            % Discrepant with expected results for now.
            
            % Subset region of interest.
            if ~isequal(size(V, 1:3), vf.span)
                V = vf.subsetField(V);
            end
            div = NaN(size(V, 1:3));
            
            switch vf.solver.diff.order
                case 1
                    % Ensure diff direction in increasing x, y, z by boolean multiplication.
                    div = squeeze(vf.diff1(V(:,:,:,1), [1 0 0], vf.solver.diff.mode)*(vf.xsp>0) + ...
                        vf.diff1(V(:,:,:,2), [0 1 0], vf.solver.diff.mode)*(vf.ysp>0) + ...
                        vf.diff1(V(:,:,:,3), [0 0 1], vf.solver.diff.mode)*(vf.zsp>0));
            end
        end
        
        function C = curl(vf, F)
            % Compute the curl of the velocity field 'F' with the accuracy
            % order set in vf.solver.diff.err_order.
            
            if size(F, 4) ~= 3
                error('Not a vector field!')
            end
            
            C = zeros(size(F));
            
            % Special case of 2D field.
            if vf.ax == 2
                nulldim = find(vf.sps == 0);
                % Right handed indices.
                vdims = mod([nulldim + 1, nulldim + 2], 3);
                vdims(vdims==0) = 3;
                % Compute curl in just the absent vector component.
                C(:,:,:,nulldim) = vf.diff(F(:,:,:,vdims(2)), vdims(1), 1) - ...
                    vf.diff(F(:,:,:,vdims(1)), vdims(2), 1);
                C(:,:,:,vdims(1)) = 0;
                C(:,:,:,vdims(2)) = 0;
                return
            end
            % 3D velocity field.
            C(:,:,:,1) = vf.diff(F(:,:,:,3), 2, 1) - ...
                vf.diff(F(:,:,:,2), 3, 1);
            C(:,:,:,2) = vf.diff(F(:,:,:,1), 3, 1) - ...
                vf.diff(F(:,:,:,3), 1, 1);
            C(:,:,:,3) = vf.diff(F(:,:,:,2), 1, 1) - ...
                vf.diff(F(:,:,:,1), 2, 1);
        end
        
        
       %%%%%%%%%%%%%%%%%%%%%%% Integral Methods %%%%%%%%%%%%%%%%%%%%%%%%
       
       % NaN entries are ignored in the cubic surface integrations.
       
       % Scalar surface integral (scalar area elements).
       function summed = intCubicSurf(vf, F)
           % Integrate on the faces of the current effective region, which
           % is a cube, the given field, vector or scalar, 'F', using a
           % scalar surface elements.
           % 
           % The indices of V are assumed to match that of X_e.
           
           if ~isequal(vf.span, size(F, 1:3))
               error('Mismatching Grids Dimensions!')
           end
           
           summed = abs(vf.ysp*vf.zsp)*(sum(F(:,vf.ascLim(1,1),:,:), 'all', 'omitnan') + ...
               sum(F(:,vf.ascLim(1,2),:,:), 'all', 'omitnan')) + ...
               abs(vf.xsp*vf.zsp)*(sum(F(vf.ascLim(2,1),:,:,:), 'all', 'omitnan') + ...
               sum(F(vf.ascLim(2,2),:,:,:), 'all', 'omitnan')) + ...
               abs(vf.xsp*vf.ysp)*(sum(F(:,:,vf.ascLim(3,1),:), 'all', 'omitnan') + ...
               sum(F(:,:,vf.ascLim(3,2),:), 'all', 'omitnan'));
           summed = squeeze(summed);
       end
       
       % Vector integral of scalar field paired with vector surface
       % elements.
       function vec = intCubicSurf_vec(vf, S)
           % Integrate on the faces of the current effective region, which
           % is a cube, the six vector surface elements multiplied by the
           % given scalar field.
           % 
           % The indices of V are assumed to match that of X_e.
           
           if ~isequal(vf.span, size(S, 1:3))
               error('Mismatching Grids Dimensions!')
           end
           
           vec = abs(vf.ysp*vf.zsp)*([-1 0 0]'*sum(S(:,vf.ascLim(1,1),:), 'all', 'omitnan') + ...
               [1 0 0]'*sum(S(:,vf.ascLim(1,2),:), 'all', 'omitnan')) + ...
               abs(vf.xsp*vf.zsp)*([0 -1 0]'*sum(S(vf.ascLim(2,1),:,:), 'all', 'omitnan') + ...
               [0 1 0]'*sum(S(vf.ascLim(2,2),:,:), 'all', 'omitnan')) + ...
               abs(vf.xsp*vf.ysp)*([0 0 -1]'*sum(S(:,:,vf.ascLim(3,1)), 'all', 'omitnan') + ...
               [0 0 1]'*sum(S(:,:,vf.ascLim(3,2)), 'all', 'omitnan'));
       end
       
       function flux = intCubicSurf_flux(vf, V, out)
           % Compute the flux on the six faces of the cube as defined by
           % the current effective region by the given vector field 'V'.
           % The indices of 'V' are assumed to match that of X_e.
           
           % output 'sum' returns the total flux over all surfaces
           % of the volume
           % output 'faces' returns the impulse flux on each face in the
           % order: [LEFT RIGHT BOTTOM TOP BACK FRONT]
           
           if ~isequal(vf.span, size(V, 1:3))
               error('Mismatching Grids Dimensions!')
           end
           
           if ~exist('out', 'var')
               out = 'sum';
           end
           
           % Left face.
           flux_left = sum(-V(:,vf.ascLim(1,1),:,1)*abs(vf.ysp)*abs(vf.zsp), 'all', 'omitnan');
           % Right face.
           flux_right = sum(V(:,vf.ascLim(1,2),:,1)*abs(vf.ysp)*abs(vf.zsp), 'all', 'omitnan');
           % Bottom face.
           flux_bottom = sum(-V(vf.ascLim(2,1),:,:,2)*abs(vf.xsp)*abs(vf.zsp), 'all', 'omitnan');
           % Top face.
           flux_top = sum(V(vf.ascLim(2,2),:,:,2)*abs(vf.xsp)*abs(vf.zsp), 'all', 'omitnan');
           % Back face.
           flux_back = sum(-V(:,:,vf.ascLim(3,1),3)*abs(vf.xsp)*abs(vf.ysp), 'all', 'omitnan');
           % Front face.
           flux_front = sum(V(:,:,vf.ascLim(3,2),3)*abs(vf.xsp)*abs(vf.ysp), 'all', 'omitnan');
           
           switch out
               case 'sum'
                    flux = flux_left + flux_right + flux_top + flux_bottom + flux_front + flux_back;
               case 'faces'
                    flux = [flux_left flux_right flux_top flux_bottom flux_front flux_back];
           end
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
           left = abs(vf.ysp)*abs(vf.zsp)*sum(crs_pd, [1 2], 'omitnan');
           % Right face.
           crs_pd = zeros(vf.span(1), vf.span(3), 3);
           crs_pd(:,:,2) = squeeze(V(:, vf.ascLim(1,2), :, 3));
           crs_pd(:,:,3) = -squeeze(V(:, vf.ascLim(1,2), :, 2));
           right = abs(vf.ysp)*abs(vf.zsp)*sum(crs_pd, [1 2], 'omitnan');
           % Bottom face.
           crs_pd = zeros(vf.span(2), vf.span(3), 3);
           crs_pd(:,:,1) = squeeze(V(vf.ascLim(2,1), :, :, 3));
           crs_pd(:,:,3) = -squeeze(V(vf.ascLim(2,1), :, :, 1));
           bottom = abs(vf.xsp)*abs(vf.zsp)*sum(crs_pd, [1 2], 'omitnan');
           % Top face.
           crs_pd = zeros(vf.span(2), vf.span(3), 3);
           crs_pd(:,:,1) = -squeeze(V(vf.ascLim(2,2), :, :, 3));
           crs_pd(:,:,3) = squeeze(V(vf.ascLim(2,2), :, :, 1));
           top = abs(vf.xsp)*abs(vf.zsp)*sum(crs_pd, [1 2], 'omitnan');
           % Back face.
           crs_pd = zeros(vf.span(1), vf.span(2), 3);
           crs_pd(:,:,1) = -squeeze(V(:, :, vf.ascLim(3,1), 2));
           crs_pd(:,:,2) = squeeze(V(:, :, vf.ascLim(3,1), 1));
           back = abs(vf.xsp)*abs(vf.ysp)*sum(crs_pd, [1 2], 'omitnan');
           % Front face.
           crs_pd = zeros(vf.span(1), vf.span(2), 3);
           crs_pd(:,:,1) = squeeze(V(:, :, vf.ascLim(3,2), 2));
           crs_pd(:,:,2) = -squeeze(V(:, :, vf.ascLim(3,2), 1));
           front = abs(vf.xsp)*abs(vf.ysp)*sum(crs_pd, [1 2], 'omitnan');
           
           crs = squeeze(left + right + bottom + top + back + front);
       end
        
    end
end