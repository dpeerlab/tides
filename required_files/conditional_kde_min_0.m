function [ estimated_points_grid, meshX, meshY, conditional_mean, xvalues] = conditional_kde_min_0(edge_data, num_points,varargin)

  

    %q = kde(transpose(edge_data),'rot');
     min_vals = min(edge_data);
     [~, y0, ~, ~, ~, x1, y1] = pairwise_visualize_matrix_min_0(edge_data, 'no_plot', varargin{:});
     %max(edge_data);
     max_vals = [max(x1), max(edge_data(:, 2))];%max(edge_data);%
     %max_vals = [max(x1), max(y1)];
     
     if(length(varargin)>0)
         for i=1:length(varargin)
             if(strcmp(varargin{i},'maxy'))
                
                 max_vals(2) = varargin{i+1};
                 
             end
             if(strcmp(varargin{i},'maxx'))
                
                 max_vals(1) = varargin{i+1};
                 
             end
             
             if(strcmp(varargin{i}, 'minx'))
                 min_vals(1) = varargin{i+1};
             end
             
             if(strcmp(varargin{i}, 'miny'))
                 min_vals(2) = varargin{i+1};
             end
             
             if (strcmp(varargin{i}, 'Limits'))
                lims_vals = varargin{i+1};
                min_vals(1) = lims_vals(1);
                min_vals(2) = lims_vals(2);
                max_vals(1) = lims_vals(3);
                max_vals(2) = lims_vals(4);
             end
             
             
         end
     end
     
% min_vals
% max_vals
     
     
     num_points = num_points + 1;
     num_points1 = 256;
     widthX = max_vals(1)/num_points1;
     widthY = max_vals(2)/num_points1;
    
     min_vals
     max_vals
     [~, estimated_points_grid, meshX, meshY] = kde2d(edge_data, 256, min_vals, max_vals);
     
     
     
%     [meshX, meshY] = meshgrid(widthX:widthX:max_vals(1), widthY:widthY:max_vals(2));  
%     meshX_skinny = reshape(meshX, 1,prod(size(meshX)));
%     meshY_skinny = reshape(meshY, 1,prod(size(meshY)));
%     estimated_points = evaluate(q,[meshX_skinny; meshY_skinny]);
%     estimated_points_grid = reshape(estimated_points, size(meshX,1), size(meshX,2));
    



    %marg_q = marginal(q, 1);
    %estimated_marginal_points = evaluate(marg_q, xvalues);
    
    
    conditional_estimated_points_grid = estimated_points_grid;
    xvalues = widthX:widthX:max_vals(1);
    yvalues = widthY:widthY:max_vals(2);
   
    conditional_mean = zeros(length(xvalues),1);
   
    
     for i = 1:length(xvalues)
       
         %conditional_estimated_points_grid(:,i) = estimated_points_grid(:,i)./estimated_marginal_points(i); 
         conditional_estimated_points_grid(:,i) = estimated_points_grid(:,i)./sum(estimated_points_grid(:,i)); 
        
         yvector = conditional_estimated_points_grid(:,i);
    
         conditional_mean(i) = dot(yvector, yvalues); 
        
     end
  
    for i = 1:length(xvalues)
        
      
       
        %estimated_points_grid(:,i) = estimated_points_grid(:,i)./max(estimated_points_grid(:,i));
        estimated_points_grid(:, i) = conditional_estimated_points_grid(:, i)./max(conditional_estimated_points_grid(:, i));
       
    end
    
  
    
%     imagesc(xvalues, yvalues, estimated_points_grid);
%     set(gca,'YDir','normal');
%     
%     figure 
%     
%     plot(xvalues, conditional_mean);
    
end

