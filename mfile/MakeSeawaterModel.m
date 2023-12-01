function velocity_distribution = MakeSeawaterModel(grid_size, dz, c0, epsilon, z0, B)
    % Define the depth vector based on the grid size and spacing
    depth = 0:dz:(grid_size(2)-1)*dz;
    
    % Initialize the velocity distribution matrix
    velocity_distribution = zeros(grid_size);
    
    % Calculate the velocity distribution based on the Munk profile
    for i = 1:grid_size(1)
        for j = 1:grid_size(2)
            z = depth(j);
            velocity_distribution(i,j) = c0 * (1 + epsilon * ((z - z0) / B)^2);
        end
    end
    velocity_distribution=velocity_distribution';
end