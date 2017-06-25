function view_2dslice_min(estimated_points_grid, meshx, meshy, meshz, labels, fix_variable, fix_value)

%estimated_points_grid(estimated_points_grid < 0.55) = 0;
%estimated_points_grid = smoothts(estimated_points_grid, 'g', 1, 0.00001);
%     for i = 2:size(estimated_points_grid, 1)-1
%         for j = 2:size(estimated_points_grid, 2)-1
%             a = estimated_points_grid(:, i-1, j);
%             b = estimated_points_grid(:, i, j);
%             c = estimated_points_grid(:, i+1, j);
%             d = estimated_points_grid(:, i, j-1);
%             e = estimated_points_grid(:, i, j+1);
%             if prod([a, b, c, d, e]) ~= 0
%                 estimated_points_grid(i, j, :) = (a + b+ c+ d  +e)/5;
%             end
%         end
%     end
    %estimated_points_grid(1:4, 1:4, 1:4)

    
    %estimated_points_grid = estimated_points_grid(estimated_points_grid > 0.6);
    if(fix_variable ==1)
        slice(meshx, meshy, meshz, estimated_points_grid, fix_value, [], []);
        showaxes('on');
        %box on
    elseif(fix_variable==2)
        slice(meshx, meshy, meshz, estimated_points_grid, [], fix_value, []);
        showaxes('on');
        %box on
    elseif(fix_variable==3)
        slice(meshx, meshy, meshz, estimated_points_grid, [], [], fix_value);
        showaxes('on');
        %box on
    end
        
    
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    set(gca, 'FontSize', 16, 'FontWeight', 'bold')
    
    xmax = max(meshx(:));
    ymax = max(meshy(:));
    zmax = max(meshz(:));
    
    
    xmin = min(meshx(:));
    ymin = min(meshy(:));
    zmin = min(meshz(:));
    
    set(gca,'XLim',[xmin xmax]);
    set(gca,'YLim',[ymin ymax]);
    set(gca,'Zlim',[zmin zmax]);
    
    
    
    
%snaps to nearest grid slice
   
%   if(fix_variable == 1)
%       fvalues = meshy(:,1,1); %anything but x here
%       fvalues = abs(fvalues - fix_value);
%       [~,index ] = min(fvalues);
%       
%       
%       slice_data = estimated_points_grid(index,:,:);
%       valuesdim1 = meshx(1,:,1);
%       valuesdim2 = meshx(1,1,:);
%       
%   elseif(fix_variable == 2)
%       fvalues = meshz(1,:,1); %anything but y here
%       fvalues = abs(fvalues-fix_value);
%       [~,index ] = min(fvalues);
%       
%       slice_data = estimated_points_grid(:,index,:);
%       valuesdim1 = meshy(:,1,1);
%       valuesdim2 = meshy(1,1,:);
%       
%   elseif(fix_variable == 3)
%       
%       fvalues = meshx(1,1,:); % anything but z here
%       fvalues = abs(fvalues - fix_value);
%       [~,index ] = min(fvalues);
%       
%       slice_data = estimated_points_grid(:,:,index);
%       valuesdim1 = meshz(:,1,1);
%       valuesdim2 = meshz(1,:,1);
%       
%   end
%   
%   
%   
%   imagesc(valuesdim1, valuesdim2, slice_data)
%   set(gca,'YDir','normal');

end

