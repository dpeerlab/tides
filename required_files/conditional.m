function [estimated_points_grid, meshX, meshY, meshZ, conditional_mean, zvalues, conditional_estimated_points_grid, estimated_points, meshX2d, meshY2d, min_stuff, stuff ] = conditional_min(edge_data, num_points,varargin)

  %q = kde(transpose(edge_data),'rot');
  %num_points = 40;
%     GridX = linspace(0,max(edge_data(:,1)),num_points);
%     GridY = linspace(0,max(edge_data(:,2)),num_points);
%     GridZ = linspace(0,max(edge_data(:,3)),num_points);
    
    %edge_data = edge_data;
    set_min = 0;
    [~, y0, ~, ~, ~, x1, y1] = pairwise_visualize_matrix_min(edge_data(:, 2:3), 'no_plot');
    

    %[~, ~, ~, ~, ~, x2, y2] = pairwise_visualize_matrix_min(edge_data(:, 2:3), 'no_plot');
    min_vals = min(edge_data);%[min(edge_data(:, 1)), min(y1), min(edge_data(:, 3))]
    %min_vals = [0, 0, 0];
    %min_vals = [min(edge_data(:, 1)), min(x1), min(y1)]
    %min_vals = [0, 2.5, 2]
    %min_vals(2) = 1;

    
    max_vals = [max(edge_data(:, 1)), max(x1), max(edge_data(:, 3))];
    %max_vals = [1, max(x1), max(y1)];
    %max_vals(1) = 0.9;
    
    %max_vals = max(edge_data);
%     min_vals
%     max_vals
%     
    
    bandwidth_yes = 0;
    if(length(varargin)>0)
         for i=1:length(varargin)
             if(strcmp(varargin{i},'maxy'))
                
                 max_vals(2) = varargin{i+1};
                 
             end
             if(strcmp(varargin{i},'maxx'))
                
                 max_vals(1) = varargin{i+1};
                 
             end
             if(strcmp(varargin{i},'maxz'))
                
                 max_vals(3) = varargin{i+1};
                 
             end
             
             
             if(strcmp(varargin{i},'miny'))
                
                 min_vals(2) = varargin{i+1};
                 
             end
             if(strcmp(varargin{i},'minx'))
                
                 min_vals(1) = varargin{i+1};
                 
             end
             if(strcmp(varargin{i},'minz'))
                
                 min_vals(3) = varargin{i+1};
                 
             end
             
             if (strcmp(varargin{i}, 'bandwidth'))
                 bandwidth_yes = 1;
                 bandwidth = varargin{i+1};
             end
             
             
             if (strcmp(varargin{i}, 'set_min'))
                 set_min =1;
             end
         end
     end
    %max_vals
    
    
    
    widthX = max_vals(1)/num_points;
    widthY = max_vals(2)/num_points;
    widthZ = max_vals(3)/num_points;
    [meshX, meshY, meshZ] = ndgrid(widthX:widthX:max_vals(1), widthY:widthY:max_vals(2), widthZ:widthZ:max_vals(3));
    %size(meshX)
    
    meshX_skinny = reshape(meshX, 1, prod(size(meshX)));
    meshY_skinny = reshape(meshY, 1, prod(size(meshY)));
    meshZ_skinny = reshape(meshZ, 1, prod(size(meshZ)));
    
    
    [meshX2d, meshY2d] = meshgrid(widthX:widthX:max_vals(1), widthY:widthY:max_vals(2));  
    meshX_skinny2d = reshape(meshX2d, 1,prod(size(meshX2d)));
    meshY_skinny2d = reshape(meshY2d, 1,prod(size(meshX2d)));
    
    
   % estimated_points = evaluate(q,[meshX_skinny; meshY_skinny; meshZ_skinny]);
   
   if bandwidth_yes == 1
       if set_min == 1
           a = kde3d(edge_data, num_points, min(edge_data), max_vals, bandwidth);           
       else
           a = kde3d(edge_data, num_points, [0, 0, 0], max_vals, bandwidth);           
       end
   else
       if set_min == 1           
           a = kde3d_min(edge_data, num_points, min_vals, max_vals);     
%            min_vals
%            max_vals
       else
           a = kde3d(edge_data, num_points, [0, 0, 0], max_vals);
       end
       %a = kde3d(edge_data, num_points, [0, 0, 0], max_vals);%, bandwidth);
   end
   
   estimated_points = a.densityEstimates();   
   xmesh = a.xmesh();
   ymesh = a.ymesh();
   zmesh = a.zmesh();
   
%    min(xmesh)
%    min(ymesh)
%    min(zmesh)
 

   [meshX, meshY, meshZ] = ndgrid(xmesh, ymesh, zmesh);
   
   estimated_points = permute(estimated_points, [2, 1, 3]);
   meshX = permute(meshX, [2, 1, 3]);
   meshY = permute(meshY, [2, 1, 3]);
   
%    estimated_points_grid = reshape(estimated_points, size(meshX,1), size(meshX,2), size(meshX,3));
    
estimated_points_grid = estimated_points;

%estimated_points_grid = evaluate(q,[meshX_skinny; meshY_skinny; meshZ_skinny]);
refinery_est = estimated_points_grid >= 0;


    %marg_q = marginal(q, [1 2]);
    %Doubly reshape?
    %estimated_marginal_points = evaluate(marg_q, [meshX_skinny2d; meshY_skinny2d]);
    %estimated_marginal_points_grid = reshape(estimated_marginal_points, size(meshX2d,1),size(meshX2d,2));
    estimated_points_grid = estimated_points_grid .* refinery_est;
    
    
    %estimated_points_grid = smooth3(estimated_points_grid);%, 'gaussian', 15);, a.bandwidth);
    

    %%%%%%%%%%%
%     std(estimated_points_grid(:))
%     refine2 = estimated_points_grid >= min(estimated_points_grid(:)) + 1/50*std(estimated_points_grid(:));
%     estimated_points_grid = estimated_points_grid .* refine2;
    %%%%%%%%%%%%%%%%%%%%%%
    
    conditional_estimated_points_grid = estimated_points_grid;% .* refinery_est;
    num_points_x = size(estimated_points_grid,1);
    
    num_points_y = size(estimated_points_grid,2); 
    
    conditional_mean = zeros(num_points_x, num_points_y);
    %zvalues = widthZ:widthZ:max_vals(3);
    zvalues = zmesh;
	%zvector = zeros(1,length(zvalues));
    
%     est_heat_new = estimated_points_grid;
%     
%     for i = 1:size(estimated_points_grid,1);
%         for j = 1:size(estimated_points_grid,2);
%             for k = 1:size(estimated_points_grid,3);
%                 if(est_heat_new(i, j, k) < 1e-10);
%                     est_heat_new(i, j, k) = 0;
%                 end
%             end
%         end
%     end
%     
%     estimated_points_grid = est_heat_new;
    
   
    
    
    
         
    %count = 0;
%     stuff = [];
%     min_stuff = [];
    for i = 1:num_points
       for j = 1: num_points 
           if sum(estimated_points_grid(i, j, :)) ~= 0 
                %sum(estimated_points_grid(i, j, :))
                conditional_estimated_points_grid(i,j,:) = estimated_points_grid(i,j,:)./sum(estimated_points_grid(i,j,:));           
                %conditional_estimated_points_grid(i,j,:) = estimated_points_grid(i,j,:)
           end
           
           
           
%            zvector = zeros(1,length(zvalues));
%            zvector(1:length(zvalues)) = conditional_estimated_points_grid(i,j,:);
            zvector = conditional_estimated_points_grid(i,j,:);
            zvector = zvector(:);
%             stuff = [stuff, max(zvector)];
%             min_stuff = [min_stuff, min(zvector)];             
            
            conditional_mean(i,j) = dot(zvector', zvalues);      
              
       end
    end
   % conditional_mean = conditional_mean * (ymesh(2) - ymesh(1)) * (xmesh(2) - xmesh(1));
    

    for i = 1:size(estimated_points_grid,1);
       for j = 1:size(estimated_points_grid,2);
           
             if max(estimated_points_grid(i, j, :)) ~= 0
                
              %estimated_points_grid(i,j,:) = abs(estimated_points_grid(i,j,:)./max(estimated_points_grid(i,j,:)));
              estimated_points_grid(i,j,:) = abs(conditional_estimated_points_grid(i,j,:)./max(conditional_estimated_points_grid(i,j,:)));                            
             end
              
       end
    end
    
    
    %estimated_points_grid = conditional_estimated_points_grid;
    
%     for i = 1:size(estimated_points_grid, 1)
%         estimated_points_grid(
%     estimated_points_grid = gauss3filter(estimated_points_grid);
    
end



% function [ estimated_points, estimated_points_grid, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d ] = conditional(edge_data, num_points,varargin)
% 
%     
% 
% %q = kde(transpose(edge_data),'rot');
% 
% 
% 
% %    num_points = 80;
% %     GridX = linspace(0,max(edge_data(:,1)),num_points);
% %     GridY = linspace(0,max(edge_data(:,2)),num_points);
% %     GridZ = linspace(0,max(edge_data(:,3)),num_points);
%     
%     
%     max_vals = max(edge_data);
%     
%     
%     a = kde3d(edge_data, num_points, [0, 0, 0], max_vals);
%     if(length(varargin)>0)
%          for i=1:length(varargin)
%              if(strcmp(varargin{i},'maxy'))
%                 
%                  max_vals(2) = varargin{i+1};
%                  
%              end
%              if(strcmp(varargin{i},'maxx'))
%                 
%                  max_vals(1) = varargin{i+1};
%                  
%              end
%              if(strcmp(varargin{i},'maxz'))
%                 
%                  max_vals(3) = varargin{i+1};
%                  
%              end
%          end
%      end
%     
%     
%     widthX = max_vals(1)/num_points;
%     widthY = max_vals(2)/num_points;
%     widthZ = max_vals(3)/num_points;
%     [meshX, meshY, meshZ] = meshgrid(widthX:widthX:max_vals(1), widthY:widthY:max_vals(2), widthZ:widthZ:max_vals(3));  
%     size(meshX)
%     
%     meshX_skinny = reshape(meshX, 1, prod(size(meshX)));
%     meshY_skinny = reshape(meshY, 1, prod(size(meshY)));
%     meshZ_skinny = reshape(meshZ, 1, prod(size(meshZ)));
%     
%     
%     [meshX2d, meshY2d] = meshgrid(widthX:widthX:max_vals(1), widthY:widthY:max_vals(2));  
%     meshX_skinny2d = reshape(meshX2d, 1,prod(size(meshX2d)));
%     meshY_skinny2d = reshape(meshY2d, 1,prod(size(meshX2d)));
%     
%     
%     %estimated_points = evaluate(q,[meshX_skinny; meshY_skinny; meshZ_skinny]);
%     estimated_points = a.densityEstimates();
%     estimates = estimated_points(:);
%     estimates(estimates < 0.000001) = 0;
%     estimated_points = reshape(estimates, size(estimated_points));
%     %estimated_points_grid1 = reshape(estimated_points, size(meshX,1), size(meshX,2), size(meshX,3));
%     estimated_points_grid = reshape(estimated_points, size(meshX,1), size(meshX,2), size(meshX,3));
%     
%     %marg_q = marginal(q, [1 2]);
%     %Doubly reshape?
%     %estimated_marginal_points = evaluate(marg_q, [meshX_skinny2d; meshY_skinny2d]);
%     %estimated_marginal_points_grid = reshape(estimated_marginal_points, size(meshX2d,1),size(meshX2d,2));
%     conditional_estimated_points_grid = estimated_points_grid;
%     num_points_x = size(estimated_points_grid,1);
%     
%     num_points_y = size(estimated_points_grid,2); 
%     
%     conditional_mean = zeros(num_points_x, num_points_y);
%     zvalues = widthZ:widthZ:max_vals(3);
%     
%     
%     
%    
%     for i = 1:size(estimated_points_grid,1);
%        for j = 1:size(estimated_points_grid,2);
%            
%    
%                 
%               estimated_points_grid(i,j,:) = estimated_points_grid(i,j,:)./max(estimated_points_grid(i,j,:));
%               
%               
%        end
%     end
%     
%     for i = 1:num_points
%        for j = 1: num_points 
% %             conditional_estimated_points_grid(i,j,:) = estimated_points_grid(i,j,:)./sum(estimated_points_grid(i,j,:)); 
% %             zvector = zeros(1,length(zvalues));
% %             zvector(1:length(zvalues)) = conditional_estimated_points_grid(i,j,:);
% %             conditional_mean(i,j) = dot(zvector, zvalues);   
% %             [value, index] = max(zvector);
% %             conditional_mean(i,j) = zvalues(index);
%               
%             
%             conditional_estimated_points_grid(i,j,:) = estimated_points_grid(i,j,:)./sum(estimated_points_grid(i,j,:));
%             %conditional_estimated_points_grid(i,j,:) = estimated_points_filtered(i,j,:)./sum(estimated_points_filtered(i,j,:));
%             %zvector = zeros(1,length(zcoords));
%             
%             estimated_points_filtered = estimated_points_grid(i,j,:)> .6;
%             estimated_points_filtered = estimated_points_filtered .* estimated_points_grid(i,j,:);
%             
%             new_weights = estimated_points_filtered./sum(estimated_points_filtered);
%             
%             %zvector(1:length(zvalues)) = new_weights;
%             zvector(1:length(zvalues)) = conditional_estimated_points_grid(i,j,:);
%             conditional_mean(i,j) = dot(zvector,zvalues); 
%             
%             
%        end
%     end
%     
% end