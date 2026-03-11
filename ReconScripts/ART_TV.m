function [recons] = ART_TV(ProblemSetup, pm_ARTTV)
    % Accurate image reconstruction from few-views and limited-angle  data in divergent-beam CT
    %       Emil Sidky 2009
    % TV Paper
    % 
    % General Variables
    A                 = ProblemSetup.A;
    projections       = ProblemSetup.projections;
    N                 = ProblemSetup.N;
    num_angles        = ProblemSetup.angle_count;
    bins              = ProblemSetup.tbins;
    hull = ProblemSetup.hull(:)';
    %method specific variables
    iterations        = pm_ARTTV.Iterations;
    GradDescSteps     = pm_ARTTV.GradDescSteps; % number of TV steps
    recon             = pm_ARTTV.initial(:);    % initialization
    alpha             = pm_ARTTV.alpha;         % dTV step size
    epsilon           = pm_ARTTV.epsilon;       % for calculating dTV
    lambda            = pm_ARTTV.lambda;
    
    recons = zeros(N*N, iterations);

    subA = cell(1,num_angles);
    subP = cell(1,num_angles);
    subV = cell(1,num_angles);
    subW = cell(1,num_angles);
    for ang = 1:num_angles
        tmpA = A((ang-1)*bins+1:ang*bins,:).*hull;
        subA{ang} = tmpA;
        subP{ang} = projections((ang-1)*bins+1:ang*bins);
        tmpV = sum(tmpA,1); tmpV = 1./tmpV; tmpV(isinf(tmpV)) = 0; tmpV = sparse(1:length(tmpV), 1:length(tmpV), tmpV, length(tmpV), length(tmpV));
        subV{ang} = tmpV;
        tmpW = sum(tmpA,2); tmpW = 1./tmpW; tmpW(isinf(tmpW)) = 0; tmpW = sparse(1:length(tmpW), 1:length(tmpW), tmpW, length(tmpW), length(tmpW)); 
        subW{ang} = tmpW;
    end

    AngleOrder        = randperm(num_angles);

    figure(1);

    for iter = 1:iterations
        prev_recon = recon;

        % (B) data projection steps
        for angle_idx = 1:num_angles
            subset = AngleOrder(angle_idx);
            A_sub = subA{subset};
            b_sub = subP{subset};
            W_sub = subW{subset};
            V_sub = subV{subset};

            recon = recon + lambda .* (V_sub * A_sub' * W_sub * (b_sub - A_sub * recon));
            recon(recon<0) = 0;
            recon(recon>3) = 3;
        end
        
        % (C)-positivity step; calculate change (to initialize distance for grad
        %           desc step
        recon(recon<0)=0;
        recon(recon>3) = 3;
        % (D)-tv gradient descent initialization
        grad_desc_size = norm(recon-prev_recon,2);
        % original exit condition
        
        % (D')grad descent
        for m = 1:GradDescSteps
            [grtv, v] = grad_TV(recon,N,epsilon);
            v = v/norm(v);
            recon = recon-alpha*grad_desc_size*v;
        end
        recons(:, iter) = recon;
imshow(full(reshape(recon,N,N)./max(recon)),'InitialMagnification','fit'); title(sprintf('ART-TV: %i of %i',iter,iterations)); pause(0.00001);
fprintf('max: %0.4f, min: %0.4f, delta: %0.4f, tv: %0.4f\n', max(recon), min(recon), norm(recon-prev_recon), grtv);
        AngleOrder        = randperm(num_angles);
    end
end

function [TV,dTV] = grad_TV(f, n, epsilon)
    f = reshape(f,n,n);
    m=n;
    %create indices in x and y directions, offset by 1 in either direction
    ind_m1 = 1:m; ind_n1 = 1:n;
    ind_m2 = [2:m 1]; ind_n2 = [2:n 1];         %assumes periodic BCs which should be ok if support is compact
    ind_m0 = [m 1:m-1]; ind_n0 = [n 1:n-1];
    diff1 = (f(ind_m2,ind_n1)-f(ind_m1,ind_n1)).^2; %difference in y-direction
    diff2 = (f(ind_m1,ind_n2)-f(ind_m1,ind_n1)).^2; %difference in x-direction
    diffttl = sqrt(diff1+diff2+epsilon^2);          %gradient norm in every pixel
    TV = sum(sum(abs(diffttl)));
    
    dTV = -1./diffttl .* (f(ind_m2,ind_n1)-2*f(ind_m1,ind_n1)+f(ind_m1,ind_n2)) + ...
            1./diffttl(ind_m0,ind_n1) .* (f(ind_m1,ind_n1)-f(ind_m0,ind_n1)) + ...
            1./diffttl(ind_m1,ind_n0) .* (f(ind_m1,ind_n1)-f(ind_m1,ind_n0));       %gradient of TV
    dTV = dTV(:);
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
function n = tv_norm(f, bins)
    % Discrete TV norm, from paper
    n = 0;
    f = reshape(f,[bins,bins]);
    for i = 2:bins
        for j = 2:bins
            n = n + sqrt((f(i,j)-f(i-1,j))^2+(f(i,j)-f(i,j-1))^2);
        end
    end
end