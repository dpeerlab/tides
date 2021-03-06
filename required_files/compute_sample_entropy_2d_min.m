%NOTE: In this code the comments are the codes that I changed and out in
%the next line of code (unless noted in specific places themselves).

function [entropy_value, valid_bit] = compute_sample_entropy_2d_min(estimated_density, num_partitions, prob_threshold)

    estimated_density_filter  = estimated_density>prob_threshold;
    estimated_density_refined = estimated_density .* estimated_density_filter;
    
    %estimated_density_refined_flat = sum(estimated_density_refined, 2);
    estimated_density_refined_flat = sum(estimated_density_refined, 1);

    [xindex] = find(estimated_density_refined_flat>0);
     
    valid_bit = 1;
    if(length(xindex)==0)
        
        entropy_value = 0;
        valid_bit = 0;
        return;
    end
    
    [numy, numx] = size(estimated_density);
    
    
    ysum_vector = zeros(numy,1);
    for i=1:length(xindex)
        
        sum_vector = zeros(numy,1);
        %sum_vector(1:numy,:) = estimated_density_refined(xindex(i),:)/sum(estimated_density_refined(xindex(i),:));
        sum_vector(1:numy,:) = estimated_density_refined(:, xindex(i))/sum(estimated_density_refined(:, xindex(i)));
        
        ysum_vector = ysum_vector + sum_vector;
        
    end

   
     y_partition_ends = 0:((numy)/num_partitions):numy;
     partition_weights = zeros(num_partitions, 1);
    
     for i=1:num_partitions
        partition_weights(i)=sum(ysum_vector(y_partition_ends(i)+1:y_partition_ends(i+1)));    
     end

     total_weight = sum(ysum_vector);
     if(total_weight == 0)
         
         entropy_value = 0;
         valid_bit = 0;
         return;
     end

     partition_weights = partition_weights./total_weight;
     
     log_partition_weights = log2(partition_weights);
     
     entropy_value = 0;
     
     for i = 1:length(partition_weights)
         if(partition_weights(i)>0)
             entropy_value = entropy_value + (partition_weights(i) * log_partition_weights(i));
         %else
             %disp('none')
         end
     end
     
     entropy_value = entropy_value* -1;
     normalization_value = log2(num_partitions);
     entropy_value = entropy_value/normalization_value';
     
end