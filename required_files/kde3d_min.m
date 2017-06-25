classdef kde3d_min
    % calculates the kernel density estimate for a 3-variate data set  
    % USING: object = kde3d(data, numBins, bandwidth, MIN_XY, MAX_XY, pointsOfEvals)
    %
    %
    % Details: 
    % Input: 1. data: the three-variate data set, must be a n by 3 matrix        
    %        2. bandwidth: a vector of three values corresponding to
    %           bandwidth in each direction is expected. The default value is
    %           [0.1, 0.1, 0.1]
    %        3. numBins: the algrotihm relies on histogram of the initial
    %                    data so the number of bins for the construction of histogram is expected
    %                    The default number of bins is 2^6.
    %                    If input provided, the number will be rounded to
    %                    the next power of 2 (as fft is used)
    %        4. MIN_XY and MAX_XY: The extremities of the grid on which the
    %                              histogram will be constructed. Default
    %                              is:
    %                              MIN_XY = min(data) - range/4
    %                              MAX_XY = max(data) + range/4
    %        5. pointsOfEvals: Points on which the density estimate will be
    %                          evaluated. 
    %                          Default is set to the grid of (numBins^3 by 3) in which each coordinate is a midpoint of some bin in the histogram 
    %                          Computation could be slow if default is set and the number of bins is large (numBins >= 2^5) .
    %                          
    %
    %%%% It is important that the pointsOfEvals lie in the range: 
    %         left = MIN_XYZ - 1/(2* (numBins - 1))* (MAX_XYZ - MIN_XYZ)
    %         right = (1- 1/(2*numBins))*(numBins /(numBins -1 )) * (MAX_XYZ - MIN_XYZ) + MIN_XYZ
    %
    %
       
    properties
        data
        numBins
        MIN_XYZ
        MAX_XYZ
        bandwidth
        tol
        pointsOfEvals
        d3dct
        d3dct_new
        effec_nBins
    end
    

    
    methods
        function obj = kde3d_min(data, numBins, MIN_XYZ, MAX_XYZ, bandwidth, tol, pointsOfEvals)
            
            % assign value to various variables if enough input not
            % provided
            if size(data, 2) ~= 3
                    fprintf('Input data must be m rows by 3 columns \n');
            else
                obj.data = data; 
            end
            
            % check of number of bins is provided
            if(nargin < 2) 
                obj.numBins = 2^7;
            else 
                obj.numBins = 2^(ceil(log2(numBins)));
            end 
            disp(obj.numBins)
            
            if(nargin == 3)
                fprintf('Max of the 2d grid not provided \n');
                return;
            end
            % check if the min and max for the grid on which histogram is
            % to be built is provided, else set it up
            if(nargin < 4)
                min_xy = min(data, [], 1);
                max_xy = max(data,[],  1);
                %ran = max_xy - min_xy;
                obj.MIN_XYZ = min_xy;% - ran/4;
                obj.MAX_XYZ = max_xy;% + ran/4;
            else
                obj.MIN_XYZ = MIN_XYZ;
                obj.MAX_XYZ = MAX_XYZ;
            end
            
            % check if bandwidth is provided
            if (nargin < 5)    
                %q = kde(transpose(obj.data),'rot');
                %obj.bandwidth = transpose(getBW(q, 1));
                sig = std(obj.data);
                obj.bandwidth = (4/5)^(1/7)*sig/(size(obj.data, 1)^(1/7));
            else
                obj.bandwidth = bandwidth;
            end
            
            
            % tolerance for truncating the effective number of coefficients
            % in the summation expression for the solution of the heat
            % equation.
            if (nargin < 6)
                obj.tol = 0;
            else
                obj.tol = tol;
            end
            
            % check if the points of evaluations is provided, else set it
            % up
            if(nargin < 7) 
                x = obj.xmesh();
                y = obj.ymesh();
                z = obj.zmesh();
                obj.pointsOfEvals = [x', y', z'];
                %[X, Y, Z] = meshgrid(x', y', z');
                %obj.pointsOfEvals = [X(:), Y(:), Z(:)];
            else
                pm = min(pointsOfEvals, [], 1); n = obj.numBins;
                pM = max(pointsOfEvals, [], 1);
                mxy = obj.MIN_XYZ;  Mxy = obj.MAX_XYZ;
                ou1 = ((pm - mxy)./(Mxy - mxy)).*(repmat((n - 1)/n,1, 3)) + repmat(1/(2*n) , 1, 3);
                ou2 = ((pM - mxy)./(Mxy - mxy)).*(repmat((n - 1)/n,1, 3)) + repmat(1/(2*n) , 1, 3);
                if ( (abs(ou1(1)) > -eps && abs(ou1(2)) > -eps && abs(ou2(1)) < 1 +eps && abs(ou2(2)) < 1 + eps))
                    obj.pointsOfEvals = pointsOfEvals;
                else
                    fprintf('points of evaluations are out of allowed range \n');
                end
            end
            
            % the 3D-DCT of the histogram stored so that it is computed
            % only once
            obj.d3dct_new = obj.botev_dct3();
            obj.d3dct = obj.dct3();
            obj.effec_nBins = obj.truncate_the_effective_number_of_bins(obj.tol);
        end
        
        % mesh in the x direction on which the histogram is constructed
        function N = xmesh(a)
            N = linspace(a.MIN_XYZ(1), a.MAX_XYZ(1), a.numBins);
        end
      
        % mesh in the y direction on which the histogram is constructed
        function N = ymesh(a)
            %N = linspace(min(a.MIN_XYZ(2), 0), a.MAX_XYZ(2), a.numBins);
            N = linspace(a.MIN_XYZ(2), a.MAX_XYZ(2), a.numBins);
        end
        
        function N = zmesh(a)
            N = linspace(a.MIN_XYZ(3), a.MAX_XYZ(3), a.numBins);
        end
        % scale xmesh to [0, 1]
        function out = scaled_xmesh(a)
           mx = a.MIN_XYZ(1); Mx = a.MAX_XYZ(1); nx = a.numBins; xmes = a.xmesh();
           out = (((xmes - mx)/(Mx - mx))*((nx-1)/(nx))) + 1/(2*nx);
        end
              
        % scale ymesh to [0, 1]
        function out = scaled_ymesh(a)
           my = a.MIN_XYZ(2); My = a.MAX_XYZ(2); ny = a.numBins; ymes = a.ymesh();
           out = (((ymes - my)/(My - my))*((ny-1)/(ny))) + 1/(2*ny);
        end
        
        % scale zmesh to [0, 1]
        function out = scaled_zmesh(a)
           mz = a.MIN_XYZ(3); Mz = a.MAX_XYZ(3); nz = a.numBins; zmes = a.zmesh();
           out = (((zmes - mz)/(Mz - mz))*((nz-1)/(nz))) + 1/(2*nz);
        end
        
        % number of rows in data
        function N = sz(a)
            N = size(a.data);
        end
        
        function r = numrows(a)
            Nbym = a.sz;
            r = Nbym(1);
        end
        
        % transform data such that data points lie in [0, 1] x [0, 1] x
        % [0,1]
        function out = transformed_data(a)
            N  = a.sz;  
            minim = repmat(a.MIN_XYZ, N(1), 1); 
            scaling = repmat(a.MAX_XYZ - a.MIN_XYZ, N(1), 1);
            int2 = (a.data - minim)./scaling;
            int3 = (a.numBins - 1)/(a.numBins);
            int4 = 1/(2*a.numBins);
            out = (int2 * int3) + int4;
        end
        
        % this is the transformation that Botev uses in his code
        function out = botev_transformed_data(a)
            N  = a.sz;  
            int0 = repmat(a.MIN_XYZ, N(1), 1); int1 = repmat(a.MAX_XYZ - a.MIN_XYZ, N(1), 1);
            out = (a.data - int0)./int1;
        end
        
        % histogram of the transformed data set
        function [h, dum] = histogram(a)
           %h = a.ndhist(a.transformed_data(), a.numBins); 
           [h, dum] = a.ndhist(a.botev_transformed_data(), a.numBins);
           %refine_h = h > 1;
           %h = h .* refine_h;
           
        end
        
        
        % histogram of botev's transformed data set
        function h = botev_histogram(a)
            h = a.ndhist(a.botev_transformed_data(), a.numBins);            
        end
        
        
        % 2D DCT of the histogram
        function d = dct3(a)
            d = a.dct3d(a.histogram());
        end
        
        % 2D DCT of Botev's histogram
        function d = botev_dct3(a)
            d = a.dct3d(a.botev_histogram());
        end
        
        % scale the bandwidth 
        function [bw] = scale_bandwidth(a)
            scaling = a.MAX_XYZ - a.MIN_XYZ;
            bw = (a.bandwidth./scaling).^2;
            %bw = bw - 0.5*bw;
        end
        
        % coefficients that go into the exponent in the summation
        % expression for the solution of the heat equation
        function [Ix] = x_expo(a)
            bw = a.scale_bandwidth(); 
            tx = bw(1);
            n = a.numBins;
            Ix = exp(-a.coeff(n)*pi^2*tx/2)';
        end
        
        % coefficients that go into the exponent in the summation
        % expression for the solution of the heat equation
        function [Ix] = x_expo2(a)
            bw = a.scale_bandwidth(); 
            tx = bw(1);
            n = a.effec_nBins; 
            Ix = exp(-a.coeff(n(1))*pi^2*tx/2)';
        end
        
        % coefficients that go into the exponent in the summation
        % expression for the solution of the heat equation
        function [Iy] = y_expo(a)
            bw = a.scale_bandwidth();
            ty = bw(2);
            n = a.numBins;
            Iy = exp(-a.coeff(n)*pi^2*ty/2)';
        end
        
        % coefficients that go into the exponent in the summation
        % expression for the solution of the heat equation
        function [Iy] = y_expo2(a)
            bw = a.scale_bandwidth();
            ty = bw(2);
            n = a.effec_nBins;
            Iy = exp(-a.coeff(n(2))*pi^2*ty/2)';
        end
        
        % coefficients that go into the exponent in the summation
        % expression for the solution of the heat equation
        function [Iz] = z_expo(a)
            bw = a.scale_bandwidth();
            tz = bw(3);
            n = a.numBins;
            Iz = exp(-a.coeff(n)*pi^2*tz/2)';
        end
        
        % coefficients that go into the exponent in the summation
        % expression for the solution of the heat equation
        function [Iz] = z_expo2(a)
            bw = a.scale_bandwidth();
            tz = bw(3);
            n = a.effec_nBins;
            Iz = exp(-a.coeff(n(3))*pi^2*tz/2)';
        end
        
        % coefficients that go into the summation
        % expression for the solution of the heat equation
        function [x] = cos_coeff_x(a)
            n = a.numBins; 
            x = (0:n-1)'*pi;
        end
        
        % coefficients that go into the summation
        % expression for the solution of the heat equation
        function [x] = cos_coeff_x2(a)
            o = a.effec_nBins; 
            x = (0:o(1)-1)'*pi;
        end
        
        % coefficients that go into the summation
        % expression for the solution of the heat equation
        function [y] = cos_coeff_y(a)
            n = a.numBins; 
            y = (0:n-1)'*pi;
        end
        
        
        % coefficients that go into the summation
        % expression for the solution of the heat equation
        function [y] = cos_coeff_y2(a)
            o = a.effec_nBins;
            y = (0:o(2)-1)'*pi;
        end
        
        % coefficients that go into the summation
        % expression for the solution of the heat equation
        function [z] = cos_coeff_z(a)
            n = a.numBins; 
            z = (0:n-1)'*pi;
        end
        
        % coefficients that go into the summation
        % expression for the solution of the heat equation
        function [z] = cos_coeff_z2(a)
            o = a.effec_nBins;
            z = (0:o(3)-1)'*pi;
        end
        
        % truncate the effective number of bins
        function [o] = truncate_the_effective_number_of_bins(a, tol)
            Ix = a.x_expo(); Iy = a.y_expo(); Iz = a.z_expo();
            if (tol == 0)
                o = [size(Ix, 1), size(Iy, 1), size(Iz, 1)];
            else
                row_x = find(Ix > 10^(-tol), 1, 'last');
                row_y = find(Iy > 10^(-tol), 1, 'last');
                row_z = find(Iz > 10^(-tol), 1, 'last');
                o = [row_x, row_y, row_z];
            end
        end
        
        % transform the points on which density is to be estimated
        function [mesh_xyz] = transform_points_of_evaluations(a)
            % pE = a.pointsOfEvals;
            N  = size(a.pointsOfEvals);  
            int0 = repmat(a.MIN_XYZ, N(1), 1); 
            int1 = repmat(a.MAX_XYZ - a.MIN_XYZ, N(1), 1);
            int2 = (a.pointsOfEvals - int0)./int1;
            int3 = (a.numBins - 1)/(a.numBins);
            mesh_xyz = (int2 * int3) + 1/(2*a.numBins);
        end
        
        % coefficients that go into the summation
        % expression for the solution of the heat equation
         function [c] = cos_matx2(a, xvec)
            x = a.cos_coeff_x2()';
            int3 = repmat(a.x_expo2()', length(xvec), 1);
            c = cos(xvec*x) .*int3; 
        end
        
        % coefficients that go into the summation
        % expression for the solution of the heat equation
        function [b] = cos_maty2(a, yvec)
            y = a.cos_coeff_y2();
            int3 = repmat(a.y_expo2(), 1, length(yvec));
            b = cos(y*yvec') .*int3; 
        end
                
        % coefficients that go into the summation
        % expression for the solution of the heat equation
        function [d] = cos_matz2(a, zvec)
            z = a.cos_coeff_z2();
            %int1 = repmat(zvec, 1, length(z));
            %int2 = repmat(z', length(zvec), 1);
            int3 = repmat(a.z_expo2()', length(zvec), 1);
            d = cos(zvec * z') .* int3;
        end
        

        % truncate the number of coefficients that go into the summation
        % expression for the solution of the heat equation
        function [bnew, row_ind_x] = truncate_x(a, xvec, tol)
            b = a.cos_matx(xvec);
            row_ind_x = find(sum(abs(b) > 10^(-tol)) > 0, 1, 'last');
            bnew = b(:, 1:row_ind_x); 
        end
        
        
        % truncate the number of coefficients that go into the summation
        % expression for the solution of the heat equation
        function [cnew, row_ind_y] = truncate_y(a, yvec, tol)
            c = a.cos_maty(yvec);
            row_ind_y = find(sum(abs(c) > 10^(-tol), 2) > 0, 1, 'last');
            cnew = c(1:row_ind_y, :); 
        end
        
        
        % truncate the number of coefficients that go into the summation
        % expression for the solution of the heat equation
        function [dnew, row_ind_z] = truncate_z(a,zvec, tol)
            d = a.cos_matz(zvec);
            row_ind_z = find(sum(abs(d) > 10^(-tol)) > 0, 1, 'last');
            dnew = d(:, 1:row_ind_z); %bnew(row_ind:end, :) = [];
        end
                      
        % get the density estimates
        function [density] = densityEstimates_NonGrid(a)
            out1 = a.transform_points_of_evaluations();
            numPoints = length(out1);
            yvec = out1(:, 2); xvec = out1(:, 1); zvec = out1(:, 3);
            bnew = a.cos_maty2(yvec); d = a.cos_matz2(zvec);
            c = a.cos_matx2(xvec);
            o = a.effec_nBins;
            dnew = a.d3dct; d1 = dnew(1:o(1), 1:o(2), 1:o(3));
            res3 = mtimesx(d1, bnew);
            c1 = repmat(c, 1, 1, o(3));
            res44 = sum(c1 .* permute(res3, [2, 1, 3]), 2);
            res55 = reshape(res44, numPoints, o(3));
            density = sum(d .* res55, 2)/prod(a.MAX_XYZ - a.MIN_XYZ);
        end

        % get the density estimates
        function [density] = densityEstimates_Grid(a)
            out1 = a.transform_points_of_evaluations();
            yvec = out1(:, 2); xvec = out1(:, 1); zvec = out1(:, 3);
            bnew = a.cos_maty2(yvec); d = a.cos_matz2(zvec); 
            c = a.cos_matx2(xvec);
            o = a.effec_nBins ;
            dnew = a.d3dct; d1 = dnew(1:o(1), 1:o(2), 1:o(3));
            res3 = mtimesx(d1, bnew);
            res4 = mtimesx(c, res3);
            res4_int = permute(res4, [1, 3, 2]);
            res5 = mtimesx(res4_int, d');
            res6 = permute(res5, [1, 3, 2]);
            density = res6/prod(a.MAX_XYZ - a.MIN_XYZ);
        end
        
        
        % reset the bandwidth if the user wants
        function a = setBandwidth(a, value)
             a.bandwidth = value;
             a.effec_nBins = a.truncate_the_effective_number_of_bins(a.tol);
        end
        
        % reset the points for evaluations if the user wants
        function a = setPointsOfEvaluations(a, value)
             a.pointsOfEvals = value;
        end
        
        % the fixed point theorem
        function [out, time] = evolve(a, t)
            N = a.numrows();
            Sum_func = a.func([0,0,2],t) + a.func([0,2,0], t) + a.func([2,0,0],t) + 2*a.func([1,1,0],t) + 2*a.func([1,0,1], t) + 2*a.func([0,1,1], t);
            time=(sqrt(64*pi*pi*pi)*(N/3)*Sum_func)^(-2/7);
            out=(t-time)/time;
            %out = time; 
        end
        
        function out=func(a, s, t)
            N = a.numrows();
            sum(s)
            if sum(s)<=4
                Sum_func=a.func([s(1)+1,s(2), s(3)],t)+a.func([s(1),s(2)+1, s(3)],t) + a.func([s(1), s(2), s(3)+1], t); 
                const=(1+1/2^(sum(s)+3/2))/3;
                time=((-1)^sum(s)*2*const*a.K(s(1))*a.K(s(2))*a.K(s(1))/N/Sum_func)^(1/(5/2+sum(s)));
                out=a.psi(s,time);
            else
                out=a.psi(s,t);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out=psi(a, s,Time)
            I = (0:a.numBins - 1).^2;
            A = a.d3dct();
            A2 = A.^2;
            % s is a vector
            w=exp(-I*pi^2*Time).*[1,.5*ones(1,length(I)-1)];
            wx=w.*(I.^s(1));
            wy=w.*(I.^s(2));
            wz=w.*(I.^s(3));
            siz = size(A2, 3);
            int = zeros(siz, 1);
            for i = 1:siz
                int(i) = wy*(A2(:, :, i))*wx';
            end
                
            out=(-1)^sum(s)*(wz*int)*pi^(2*sum(s));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function t = calc_time(a)
            t=fzero( @(t)(t-a.evolve(t)),[0,0.1]);
        end
        
        % reset the tolerance level for truncation if desired
        function a = setTolerance(a, value)
             a.tol = value; 
             a.effec_nBins = a.truncate_the_effective_number_of_bins(a.tol);
        end
        
        % get the evaluation points 
         function x = getXCoordinates(a, str1)
             coords = a.pointsOfEvals;
             if (strcmp(str1, 'nongrid'))
                x = coords(:, 1);
             elseif (strcmp(str1, 'grid') || strcmp(str1, 'default'))
                [c1, ~, ~] = meshgrid(coords(:, 1), coords(:, 2), coords(:, 3)); 
                x = permute(c1, [2, 1, 3]);
             end
         end
         
         function y = getYCoordinates(a, str1)
             coords = a.pointsOfEvals;
             if (strcmp(str1, 'nongrid'))
                y = coords(:, 2);
             elseif (strcmp(str1, 'grid') || strcmp(str1, 'default'))
                [~, c2, ~] = meshgrid(coords(:, 1), coords(:, 2), coords(:, 3)); 
                y = permute(c2, [2, 1, 3]);
             end
         end
         
         function z = getZCoordinates(a, str1)
             coords = a.pointsOfEvals;
             if (strcmp(str1, 'nongrid'))
                z = coords(:, 3);
             elseif (strcmp(str1, 'grid') || strcmp(str1, 'default'))
                [~, ~, z] = meshgrid(coords(:, 1), coords(:, 2), coords(:, 3)); 
             end
         end
        
        
        
        % process the DCT output by multiplying to it the exponential
        % portion
        function o = processedDCT(a)
            bw = a.scale_bandwidth();
            t_x = bw(1); 
            t_y = bw(2); 
            t_z = bw(3);
            xy_part = exp(-(0:a.numBins-1)'.^2*pi^2*t_x/2)*exp(-(0:a.numBins-1).^2*pi^2*t_y/2);
            z_part1 = reshape(exp(-(0:a.numBins-1).^2*pi^2*t_z/2), 1, 1, a.numBins);
            o1 = mtimesx(xy_part, z_part1);
            o = o1 .* a.d3dct; 
        end
        
        % compute the density estimates at the mid point of the histogram
        function [density] = densityEstimates_Fast(a)
            out1 = a.processedDCT(); 
            scaling = a.MAX_XYZ - a.MIN_XYZ;
            density=a.idct3d(out1)*(numel(out1)/prod(scaling));
        end
                
        
        % compute the density on the basis of what the user wants:
        function[density] = densityEstimates(a, varargin)
             
             if (isempty(varargin))
                 density = a.densityEstimates_Fast();                
             else
                if (strcmp(varargin{1}, 'nongrid'))
                    density = a.densityEstimates_NonGrid();
                elseif strcmp(varargin{1}, 'grid')
                    density = a.densityEstimates_Grid();
                else
                    fprintf('please specify grid vs nongrid option for the points of evaluations');
                    return;
                end
            end     
        end
    end
    
    
    
    %%%%%%%%%%%% SOME STATIC FUNCTIONS REQUIRED %%%%%%%%%%
    methods(Static = true)        
        % coefficients required in the summation representation of heat eq.
        % solution
        function I = coeff(n)
            I=(0:n-1).^2;
        end
        
        function out=K(s)
            out=(-1)^s*prod((1:2:2*s-1))/sqrt(2*pi);
        end
        
        function [out] = s(n)
            out = n;
        end
        
        function [binned_data, dum]=ndhist(data,M)
            % this function computes the histogram
            % of an n-dimensional data set;
            % 'data' is nrows by n columns
            % M is the number of bins used in each dimension
            % so that 'binned_data' is a hypercube with
            % size length equal to M;
            [nrows,ncols]=size(data);
            bins=zeros(nrows,ncols);
            for i=1:ncols
                [dum,bins(:,i)] = histc(data(:,i),(0:1/M:1),1);
                bins(:,i) = min(bins(:,i),M);
            end
            % Combine the  vectors of 1D bin counts into a grid of nD bin
            % counts.
            binned_data = accumarray(bins(all(bins>0,2),:),1/nrows,M(ones(1,ncols)));
            %binned_data = accumarray(bins(all(bins>10,2),:),1/nrows,M(ones(1,ncols)));
        end
        
        function data=dct3d(data)
            % computes the 3 dimensional discrete cosine transform of data
            % data is an nd cube
            [nrows,ncols, nheights]= size(data);
            if nrows~=ncols || nrows ~= nheights || ncols ~= nheights
                error('data is not a square array!')
            end
            % Compute weights to multiply DFT coefficients
            w = [1;2*(exp(-1i*(1:nrows-1)*pi/(2*nrows))).'];
            weight= w(:,ones(1,ncols), ones(1, nheights));
            int0 = dct1d(data);
            int1 = dct1d(permute(int0, [2, 3, 1]));
            int2 = dct1d(permute(int1, [2, 3, 1]));
            data = permute(int2, [2, 3, 1]);
            function transform1d=dct1d(x)
                % Re-order the elements of the columns of x
                x(:, :, :) = [x(1:2:end, :, :); x(end:-2:1, :, :)];
       
                % Multiply FFT by weights:
                transform1d = real(weight.* fft(x));
            end
        end
        
        
        function data = idct3d(data)
            % computes the 2 dimensional inverse discrete cosine transform
            [nrows,ncols, nheights]=size(data);
            % Compute weights
            w = exp(1i*(0:nrows-1)*pi/(2*nrows)).';
            weights=w(:,ones(1,ncols), ones(1, nheights));
            int0 = idct1d(data);
            int1 = idct1d(permute(int0, [2, 3, 1]));
            int2 = idct1d(permute(int1, [2, 3, 1]));
            data = permute(int2, [2, 3, 1]);
            %data=idct1d(idct1d(data)');
            function out=idct1d(x)
                y = real(ifft(weights.*x));
                out = zeros(nrows,ncols, nheights);
                out(1:2:nrows,:, :) = y(1:nrows/2,:, :);
                out(2:2:nrows,:, :) = y(nrows:-1:nrows/2+1,:, :);
                %out(:, :, :) = [x(1:2:end, :, :); x(end:-2:1, :, :)];
            end
        end

        
        function data=new_dct3d(data)
            % computes the 3 dimensional discrete cosine transform of data
            % data is an nd cube
            [nrows,ncols, nheights]= size(data);
            if nrows~=ncols || nrows ~= nheights || ncols ~= nheights
                error('data is not a square array!')
            end
            % Compute weights to multiply DFT coefficients
            w = [1;2*(exp(-1i*(1:nrows-1)*pi/(2*nrows))).'];
            weight= w(:,ones(1,ncols), ones(1, nheights));
            int0 = dct1d(data);
            int1 = dct1d(permute(int0, [2, 3, 1]));
            int2 = dct1d(permute(int1, [2, 3, 1]));
            data = permute(int2, [2, 3, 1]);
            function transform1d=dct1d(x)
                % Re-order the elements of the columns of x
                x(:, :, :) = [x(1:2:end, :, :); x(end:-2:1, :, :)];     
                % Multiply FFT by weights:
                transform1d = real(weight.* fft(x));
            end
        end
    end
end