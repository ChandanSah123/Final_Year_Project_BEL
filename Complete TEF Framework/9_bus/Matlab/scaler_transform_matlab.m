function scaled_data = scaler_transform_matlab(raw_data, data_min, data_max, feature_range)
    % Implements the MinMaxScaler forward transformation in MATLAB.
    % raw_data: Input data (e.g., [N x M], where N is samples, M is features)
    % data_min: Minimum values used during fitting (e.g., [1 x M])
    % data_max: Maximum values used during fitting (e.g., [1 x M])
    % feature_range: [min_range, max_range] (e.g., [-1, 1])
    
    range_min = feature_range(1);
    range_max = feature_range(2);
    
    % Ensure data_min and data_max are row vectors for element-wise operation
    data_min = data_min(:)'; 
    data_max = data_max(:)';
    
    scaled_data = (raw_data - data_min) ./ (data_max - data_min) * (range_max - range_min) + range_min;
end