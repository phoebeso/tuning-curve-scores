function hd_score = calculate_hd_score(direction,spiketrain,dt,n_hd_bins)
    dirVec = 0:2*pi/n_hd_bins:2*pi-2*pi/n_hd_bins;
    hd_count = zeros(n_hd_bins,2); % col 1 = num of times hd bin accessed, col 2 = num of spikes total

    for i = 1:numel(direction)
    
        % figure out the hd index
        [~, idx] = min(abs(direction(i)-dirVec));
        hd_count(idx,1) = hd_count(idx,1) + 1;
        hd_count(idx,2) = hd_count(idx,2) + spiketrain(i);
  
    end
    
    hd_rates = hd_count(:,2)./(hd_count(:,1).*dt); % vector with spikerates for each hd bin 
    smooth = hamming(10);
    smooth_hd = conv(hd_rates,smooth,'same');
    
    mean_vector_length = norm(smooth_hd)/n_hd_bins; % calculates the mean vector length = hd score
    hd_score = mean_vector_length; % hd score