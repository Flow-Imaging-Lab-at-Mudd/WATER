classdef CylVF < handle
    % Cylindrical parameterization of a velocity field.
    
    % Imported fields.
    properties
        
        % 4D matrix of position vectors on grid.
        X
        % 4D matrix of velocity vectors on grid.
        U
        
        % Grid dimension. Note that the indices on grid, due to meshgrid
        % convension, is given as (y, x, z), where x represent the number
        % of distinct values of x, etc.. Row vector.
        dims
        
        % Lower and upper bounds of x values.
        rbounds
        % Uniform spacing in x.
        rsp
        % Lower and upper bounds of y values.
        tbounds
        % Uniform spacing y.
        tsp
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
        
        methods (Static)
            function vf = importCmps(r, theta, z, ur, ut, uz, minimal)
                % Import velocity field by components of position and velocity,
                % as 3D arrays.
                
                dims = size(r);
                % Pack positions and velocities compactly as 3-vectors in extra dimension.
                X = zeros(dims);
                X(:,:,:,1) = r;
                X(:,:,:,2) = theta;
                X(:,:,:,3) = z;
                U = zeros(dims);
                U(:,:,:,1) = ur;
                U(:,:,:,2) = ut;
                U(:,:,:,3) = uz;
                
                vf = CylVF(X, U, minimal);
            end
        end
        
        methods
            function vf = CylVF(X, U, minimal)
                
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
                
                vf.rbounds = [X(1,1,1,1) X(1,end,1,1)];
                % Assume data is uniformly spaced in position.
                % Allowing for lower dimensional grid by defaulting sp to 0.
                vf.rsp = 0;
                if dims(2) > 1
                    vf.rsp = X(1,2,1,1) - X(1,1,1,1);
                end
                vf.tbounds = [X(1,1,1,2) X(end,1,1,2)];
                vf.tsp = 0;
                if dims(1) > 1
                    vf.tsp = X(2,1,1,2) - X(1,1,1,2);
                end
                vf.zbounds = [X(1,1,1,3) X(1,1,end,3)];
                vf.zsp = 0;
                if dims(3) > 1
                    vf.zsp = X(1,1,2,3) - X(1,1,1,3);
                end
                vf.sps = [vf.rsp vf.tsp vf.zsp];
                vf.bounds = [vf.rbounds; vf.tbounds; vf.zbounds];
                vf.lims = vf.bounds;
                
                vf.range = [ones(3, 1) vf.dims'];
                vf.span = (vf.range*[-1; 1])' + 1;
                
                % Initialize innate quantities and customized settings.
                vf.initPropertyStructs()
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
            end
            
            function K = kineticEnergy(vf)
                % Compute the kinetic energy of the effective region with
                % definitional formula. NaN values of velocity are ignored in
                % sum.
                
                K = 1/2*vf.fluid.density * ...
                    sum(vf.X(:,:,:,1)*vf.tsp*vf.rsp*vf.zsp * vf.scale.len^5 .* ...
                    sum(vf.U.^2, 4), 'all', 'omitnan');
            end
            
        end
end