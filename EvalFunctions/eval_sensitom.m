function overall_eval = eval_sensitom(recons, N, iterations, EvalParameters)
    % major targets
    phantom_insert_locations = EvalParameters.Insert_Locations;
    inserts = fieldnames(phantom_insert_locations);
    insert_truth = EvalParameters.Insert_Truth;
    phantom_insert_radius = EvalParameters.InsertRadius;

    % build template fields
    template = struct();
    for ph_idx = 1:length(inserts)
        mat = inserts{ph_idx};
        template.(mat + "_truth") = [];
        template.(mat + "_recon") = [];
        template.(mat + "_error") = [];
    end
    
    overall_eval = repmat(template, length(iterations), 1);

    for iter = 1:iterations
        recon = reshape(recons(:,iter),N,N);
        % figure(); imshow(recon./max(recon(:)), 'InitialMagnification','fit');  

% hold on; 
        eval = struct();
        for ph_idx = 1:length(inserts)
            insert_material = inserts{ph_idx};
            insert_loc = phantom_insert_locations.(insert_material);
    
            % draw_circle(insert_loc(1), insert_loc(2), phantom_insert_radius,insert_material);
            intensity = circle_intensity(recon, insert_loc(1), insert_loc(2), phantom_insert_radius);
            
            error = abs(insert_truth(insert_material) - intensity) / (insert_truth(insert_material) + 1e-6);
            eval.(strcat(insert_material,"_truth")) = insert_truth(insert_material);
            eval.(strcat(insert_material,"_recon")) = intensity;
            eval.(strcat(insert_material,"_error")) = error;
        end
        overall_eval(iter) = eval;

    end
    
end

function draw_circle(xc, yc, r,label)
    theta = linspace(0, 2*pi, 200);
    x = xc + r*cos(theta);
    y = yc + r*sin(theta);
    plot(x, y, 'Color', 'r', 'LineWidth', .5);
    text(xc+r, yc+r, label);
end

function avg_intensity = circle_intensity(A, xc, yc, r)
    [x,y] = meshgrid(1:size(A,1), 1:size(A,1));
    mask = (x - xc).^2 + (y - yc).^2  <= r^2;
    avg_intensity = mean(A(mask));
end