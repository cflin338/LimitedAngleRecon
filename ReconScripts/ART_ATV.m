function recons = ART_ATV(ProblemSetup, pm_ARTATV)
    % A limited-angle CT reconstruction method based on anisotropic TV minimization
    %       Zhiqiang Chen 2013
    %       SART + ATV

    % General Variables
    PseuTarget = ProblemSetup.PseuTarget;
    A                 = ProblemSetup.A;
    projections       = ProblemSetup.projections;
    N                 = ProblemSetup.N;
    thetas            = ProblemSetup.angles; %angles
    bins              = ProblemSetup.tbins;
    num_angles        = ProblemSetup.angle_count;
        formatted_proj = reshape(projections, num_angles, bins);
        sim_b = ProblemSetup.projections_sim;
        sim_b_formatted = reshape(sim_b, num_angles,bins);
    % method specific variables
    recon             = pm_ARTATV.initial;
    N_tv              = pm_ARTATV.N_tv;
    alphas            = pm_ARTATV.alphas; %atv directions
    iterations        = pm_ARTATV.iterations;
    lambda            = pm_ARTATV.lambda;
    tvStep = pm_ARTATV.tvStep;
    subA = cell(1,num_angles);
    subP = cell(1,num_angles);
    subV = cell(1,num_angles);
    subW = cell(1,num_angles);

    sim_b = ProblemSetup.projections_sim;
    subP_sim = cell(1,num_angles);
    % figure(2);
    for ang = 1:num_angles
        tmpA = A((ang-1)*bins+1:ang*bins,:);
        % for binID = 1:bins
        %     imagesc(reshape(tmpA(binID,:),N,N)); pause(0.001);
        % end
        subA{ang} = tmpA;
        subP{ang} = projections((ang-1)*bins+1:ang*bins);
        subP_sim{ang} = sim_b((ang-1)*bins+1:ang*bins);
        tmpV = sum(tmpA,1); tmpV = 1./tmpV; tmpV(isinf(tmpV)) = 0; tmpV = sparse(1:length(tmpV), 1:length(tmpV), tmpV, length(tmpV), length(tmpV));
        subV{ang} = tmpV;
        tmpW = sum(tmpA,2); tmpW = 1./tmpW; tmpW(isinf(tmpW)) = 0; tmpW = sparse(1:length(tmpW), 1:length(tmpW), tmpW, length(tmpW), length(tmpW)); 
        subW{ang} = tmpW;
        % imagesc(reshape(sum(tmpA,1),N,N)); title(ang); pause(0.0001); 
    end
    
    recons = zeros(N*N, iterations);

    AngleOrder        = randperm(num_angles);
    figure(1);
    for iter = 1:iterations
        prev_recon = recon;

        for angle_idx = 1:num_angles
            subset = AngleOrder(angle_idx);
            A_sub = subA{subset};
            b_sub = subP{subset};
            W_sub = subW{subset};
            V_sub = subV{subset};

            recon = recon + lambda .* (V_sub * A_sub' * W_sub * (b_sub - A_sub * recon));
            recon(recon<0) = 0;
            
            % imagesc(full(reshape(recon,N,N))); title(angle_idx); pause(0.0001);
        end
        
        % ATV gradient descent
        d = norm(recon-prev_recon,2);
         
        for n = 1:N_tv
            % gradient of ATV cost function
            [~,dtv] = grad_ATV(recon, thetas, alphas); %this is humphry's implementation
            dtv = dtv / norm(dtv,2);
            recon = recon - tvStep * d * dtv;
            % imshow(full(reshape(recon,N,N)), 'InitialMagnification','fit'); title(angle_idx); pause(0.0001);
        end

        % ------------
% epsilon = 0.1; grad_desc_size = d; alpha = 0.6;
% for n = 1:N_tv
% v = dTV(recon,N,epsilon);
% v = v/norm(v);
% recon = recon-alpha*grad_desc_size*v;
% % imagesc(full(reshape(recon,N,N))); title(n); pause(0.0001);
% end
        % ------------
        % shuffle angle order
        AngleOrder        = randperm(num_angles);
        % store reconstruction
        recons(:,iter) = recon';
        imshow(full(reshape(recon,N,N)./max(recon)), 'InitialMagnification','fit'); title(iter); pause(0.0001);
        disp([min(recon), max(recon), norm(PseuTarget - recon)])
    end

    
end



function w = ChenWeights(thetas,alphas)
%--------------------------------------------------------------------------
%       thetas: proj angles
%       alphas: directions for ATV, independent of
%                      thetas; say [0,45,90,135]
%--------------------------------------------------------------------------
% verified via TDHumphries:
%       https://github.com/TDHumphries/PSARTSUP/blob/master/core/reconstruction/compute_omega.m
% simplified version, only calculates for the very center pixel
% assums parallel beam, rotation around center pixel
% assume angles are in degrees, not radians
w = zeros(numel(alphas),1);
for alpha_idx = 1:numel(alphas)
    alpha = alphas(alpha_idx);
    w_i = 0;
    for theta = thetas
        dtheta = abs(alpha-theta);
        delta = min(dtheta, 180-dtheta);
        k = cos((90-delta)/180*pi);
        w_i = w_i + k;
    end
    w(alpha_idx,1) = w_i;
end

% normalize w_i's; guarantees ATV minimization has equal strength to TV
% minimization
w = w/norm(w,2);
end


function [TV,dTV] = grad_ATV(f, thetas, alphas)
% https://github.com/TDHumphries/PSARTSUP/blob/6bf966ae29561322baacbaa36de7da95f0646832/core/reconstruction/grad_ATV.m
    epsilon = 1e-6;

    [m, n] = size(f);
    TV = 0;
    dTV = zeros(size(f));
    omegas = zeros(1, numel(alphas));
    
    %computation of weighting parameters (see Section 2.4 of paper)
    for alpha_i = 1:numel(alphas)
        alpha = alphas(alpha_i);
        omegas(alpha_i) = compute_omega(thetas, alpha);
    end
    % See equations 17, 18
    omegas = omegas/norm(omegas, 2);
    
    %for each direction
    for alpha_i = 1:numel(alphas)
        alpha = alphas(alpha_i);
        omega = omegas(alpha_i);
        
        % Given alpha, calculate the unit vector e_alpha
        e_alpha = [cos(alpha/180*pi), sin(alpha/180*pi)]; 
        e_alpha_x = e_alpha(1); e_alpha_y = e_alpha(2); 
        
        ind_m0 = [m 1:m-1]; ind_n0 = [n 1:n-1];
        ind_m2 = [2:m 1]; ind_n2 = [2:n 1]; %assumes periodic BCs which should be ok if support is compact        
        ind_m1 = 1:m; ind_n1 = 1:n;
        
        % m,n
        diffx_1 = f(ind_m2, ind_n1) - f(ind_m1, ind_n1);
        diffy_1 = f(ind_m1, ind_n2) - f(ind_m1, ind_n1);
        
        % m-1, n
        diffx_0 = f(ind_m1, ind_n1) - f(ind_m0, ind_n1);
        diffy_0 = f(ind_m0, ind_n2) - f(ind_m0, ind_n1);        
        
        % m, n-1
        diffx_2 = f(ind_m2, ind_n0) - f(ind_m1, ind_n0);
        diffy_2 = f(ind_m1, ind_n1) - f(ind_m1, ind_n0);
        
        %dot product
        diff_1 = e_alpha_x*diffx_1 + e_alpha_y*diffy_1 + epsilon;
        diff_0 = e_alpha_x*diffx_0 + e_alpha_y*diffy_0 + epsilon;
        diff_2 = e_alpha_x*diffx_2 + e_alpha_y*diffy_2 + epsilon;        
        
        %norm of the above
        diffttl_1 = abs(diff_1);
        diffttl_0 = abs(diff_0);
        diffttl_2 = abs(diff_2);
        
        
        TV = TV + (omega * sum(sum(diffttl_1)));
        
        %gradient calculation
        alpha_i_dTV = ( (1./diffttl_1) .* (diff_1) * (-e_alpha_x - e_alpha_y) ) + ...
                      ( (1./diffttl_0) .* (diff_0) * e_alpha_x ) + ...
                      ( (1./diffttl_2) .* (diff_2) * e_alpha_y );


        dTV = dTV + (omega * alpha_i_dTV);        
    end    
end

%% unused

function v = dATV_chen(f, n, alphas, weights,epsilon)
%--------------------------------------------------------------------------
% f the current estimate, n size of original img
% alphas: directions used for ATV
% weights: weights to be used, 1-1 corresponding with alphas
% epsilon: calculating dtv to avoid div0 
%--------------------------------------------------------------------------
    % calculated with TDHumphy's implementation: 
    %https://github.com/TDHumphries/PSARTSUP/blob/master/core/reconstruction/grad_ATV.m
    f = reshape(f, [n,n]);
    v = zeros(n,n);
    for i = length(alphas)
        w = weights(i);
        alpha = alphas(i)/180*pi;
        e_x = cos(alpha/180*pi);
        e_y = sin(alpha/180*pi);
        
        ind_m0 = [n 1:n-1]; ind_n0 = [n 1:n-1];
        ind_m2 = [2:n 1]; ind_n2 = [2:n 1]; %assumes periodic BCs which should be ok if support is compact        
        ind_m1 = 1:n; ind_n1 = 1:n;
        
        % m,n
        diffx_1 = f(ind_m2, ind_n1) - f(ind_m1, ind_n1);
        diffy_1 = f(ind_m1, ind_n2) - f(ind_m1, ind_n1);
        
        % m-1, n
        diffx_0 = f(ind_m1, ind_n1) - f(ind_m0, ind_n1);
        diffy_0 = f(ind_m0, ind_n2) - f(ind_m0, ind_n1);        
        
        % m, n-1
        diffx_2 = f(ind_m2, ind_n0) - f(ind_m1, ind_n0);
        diffy_2 = f(ind_m1, ind_n1) - f(ind_m1, ind_n0);
        
        %dot product
        diff_1 = e_x*diffx_1 + e_y*diffy_1 + epsilon;
        diff_0 = e_x*diffx_0 + e_y*diffy_0 + epsilon;
        diff_2 = e_x*diffx_2 + e_y*diffy_2 + epsilon;        
        
        %norm of the above
        diffttl_1 = abs(diff_1);
        diffttl_0 = abs(diff_0);
        diffttl_2 = abs(diff_2);
        
        alpha_i_dTV = ( (1./diffttl_1) .* (diff_1) * (-e_x - e_y) ) + ...
                      ( (1./diffttl_0) .* (diff_0) * e_x ) + ...
                      ( (1./diffttl_2) .* (diff_2) * e_y );
                  
        v = v + (w * alpha_i_dTV);  
        
    end
    v = v(:);
end

function v = dATV_Jin(f, n, angles, weights, epsilon)
% this is from this paper: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=5874180
% paper picks fixed values of A, B (=1)
% implementation has varying A, B based on angle
    f = reshape(f, [n,n]);
    if nargin==4
        epsilon=1e-8;
    end
    %or fixed A values
    v = zeros(n,n);
%     for theta = angles
    for i = length(angles)
        w = weights(i);
        % theta = angles(i);
        % A = max(cos(theta/180*pi),.001);
        A = 1;
        % B = max(sin(theta/180*pi),.001);
        B = .1;
        for s = 2:size(f,1)-1
            for t = 2:size(f,2)-1
                v(s,t) = w * ((A*(f(s,t)-f(s-1,t))+B*(f(s,t)-f(s,t-1))) / ...
                        sqrt(epsilon + A*(f(s,t)-f(s-1,t))^2+B*(f(s,t)-f(s,t-1))^2) - ...
                        A*(f(s+1,t)-f(s,t))/sqrt(epsilon+A*(f(s+1,t)-f(s,t))^2+B*(f(s+1,t)-f(s+1,t-1))^2) - ...
                        B*(f(s,t+1)-f(s,t))/sqrt(epsilon+A*(f(s,t+1)-f(s,t))^2+B*(f(s,t+1)-f(s-1,t+1))^2));
            end
        end
        
    end
    v = v(:);
end

function [omega] = compute_omega(thetas, alpha)
%--------------------------------------------------------------------------
%       thetas: proj angles
%       alphas: directions for ATV, independent of
%                      thetas; say [0,45,90,135]
%--------------------------------------------------------------------------
% verified via TDHumphries:
%       https://github.com/TDHumphries/PSARTSUP/blob/master/core/reconstruction/compute_omega.m
    omega = 0;
    
    for theta_j = 1:numel(thetas)
        theta = thetas(theta_j);
  
        theta = mod(theta,360); alpha = mod(alpha,360);
        if theta > 180
            theta = theta - 180;
        end
        if alpha > 180
            alpha = alpha - 180;
        end

        br_theta = theta + 90;% changed this for MIRT angular convention
        if (br_theta > 180)
            br_theta = br_theta - 180;
        end
        
        omega = omega + abs(sin(alpha*pi/180 - br_theta*pi/180));
    end    
end



function v = dTV(f, n, epsilon)
    % gradient of TV cost function, from paper
    f = reshape(f, [n,n]);
    v = zeros(n,n);
    for s = 2:size(f,1)-1
        for t = 2:size(f,2)-1
            v(s,t) = (2*(f(s,t)-f(s-1,t))+2*(f(s,t)-f(s,t-1))) / ...
                    sqrt(epsilon + (f(s,t)-f(s-1,t))^2+(f(s,t)-f(s,t-1))^2) - ...
                    2*(f(s+1,t)-f(s,t))/sqrt(epsilon+(f(s+1,t)-f(s,t))^2+(f(s+1,t)-f(s+1,t-1))^2) - ...
                    2*(f(s,t+1)-f(s,t))/sqrt(epsilon+(f(s,t+1)-f(s,t))^2+(f(s,t+1)-f(s-1,t+1))^2);
        end
    end
    v = v(:);
end