% function [entropy_value, valid_bit] = compute_sample_entropy(estimated_density, num_partitions, prob_threshold)
% 
%     %refine density by taking out values from outliers. 
%     estimated_density_filter  = estimated_density>prob_threshold;
%     estimated_density_refined = estimated_density .* estimated_density_filter;
%     
%  
%     
%     estimated_density_refined_flat = sum(estimated_density_refined, 3);
%     [xindex, yindex] = find(estimated_density_refined_flat>0);
%     valid_bit = 1;
%     
%        
%    
%     %if none are found just return 0
%     if(length(xindex)==0)
%         entropy_value = 0;
%         valid_bit = 0;
%         return;
%     end
%     
%     %size of the estimated density
%     [numx, numy, numz] = size(estimated_density);
%     
%     %number of x y and z values 
%     zsum_vector = zeros(numz,1);
%     
%    
%     for i=1:length(xindex)
%         
%         sum_vector = zeros(numz,1);
%         sum_vector(1:numz,:) = estimated_density_refined(xindex(i),yindex(i),:);
%         zsum_vector = zsum_vector + sum_vector;
%         
%     end
% 
%    %partition z axis into the correct number of partitions
%      z_partition_ends = 0:((numz)/num_partitions):numz;
%      partition_weights = zeros(num_partitions, 1);
%      
%     %what proportion of values are in each partition
%      for i=1:num_partitions
%             
%         partition_weights(i)=sum(zsum_vector(z_partition_ends(i)+1:z_partition_ends(i+1)));
%         
%      end
% 
%      %total weight
%      total_weight = sum(zsum_vector);
%      if(total_weight == 0)
%          entropy_value = 0;
%          valid_bit = 0;
%          return;
%      end
%      partition_weights = partition_weights./total_weight;
%      
%      %formula for entropy of that vector
%      log_partition_weights = log2(partition_weights);
%      
%      entropy_value = 0; 
%      
%      for i = 1:length(partition_weights)
%          if(partition_weights(i)>0)
%              entropy_value = entropy_value + (partition_weights(i) * log_partition_weights(i));
%          end
%      end
%      
%      entropy_value = entropy_value* -1;
%      normalization_value = log2(num_partitions);
%      entropy_value = entropy_value/normalization_value;
%      
% end
% 
function [entropy_value, valid_bit] = compute_sample_entropy(estimated_density, num_partitions, prob_threshold)

    %refine density by taking out values from outliers. 
    estimated_density_filter  = estimated_density>prob_threshold;
    estimated_density_refined = estimated_density .* estimated_density_filter;
    
 
    
    estimated_density_refined_flat = sum(estimated_density_refined, 3);
    [yindex, xindex] = find(estimated_density_refined_flat>0);
    valid_bit = 1;
    
    %if none are found just return 0
    if(length(xindex)==0)
        entropy_value = 0;
        valid_bit = 0;
        return;
    end
    
    %size of the estimated density
    [numy, numx, numz] = size(estimated_density);
    
    %number of x y and z values 
    zsum_vector = zeros(numz,1);    
   
    for i=1:length(xindex)
        sum_vector = zeros(numz,1);
        %THIS IS THE LINE I CHANGED... THE CONTRIBUTION FROM EACH xindex
        %HAS TO BE EQUAL
        sum_vector(1:numz,:) = estimated_density_refined(yindex(i),xindex(i),:)/sum(estimated_density_refined(yindex(i),xindex(i),:));
        zsum_vector = zsum_vector + sum_vector;        
    end

     %partition z axis into the correct number of partitions
     z_partition_ends = 0:((numz)/num_partitions):numz;
     partition_weights = zeros(num_partitions, 1);
     
     %what proportion of values are in each partition
     for i=1:num_partitions            
        partition_weights(i)=sum(zsum_vector(z_partition_ends(i)+1:z_partition_ends(i+1)));        
     end

     %total weight
     total_weight = sum(zsum_vector);
     if(total_weight == 0)
         entropy_value = 0;
         valid_bit = 0;
         return;
     end
     partition_weights = partition_weights./total_weight;
     
     %formula for entropy of that vector
     log_partition_weights = log2(partition_weights);
     
     entropy_value = 0; 
     
     for i = 1:length(partition_weights)
         if(partition_weights(i)>0)
             entropy_value = entropy_value + (partition_weights(i) * log_partition_weights(i));
         end
     end
     
     entropy_value = entropy_value* -1;
     normalization_value = log2(num_partitions);
     entropy_value = entropy_value/normalization_value;     
end