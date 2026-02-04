function [recon,metrics, samples] = ART_TV(ProblemSetup, pm_ARTTV)
    % Accurate image reconstruction from few-views and limited-angle  data in divergent-beam CT
    %       Emil Sidky 2009
    % TV Paper
    % 
    % General Variables
    A                 = ProblemSetup.A;
    projections       = ProblemSetup.projections;
    img               = ProblemSetup.img;
    N                 = ProblemSetup.N;
    proj_order = ProblemSetup.proj_order;
    
    %method specific variables
    iterations        = pm_ARTTV.Iterations;
    NData             = pm_ARTTV.NData;
    exit_error        = pm_ARTTV.exit_error;
    GradDescSteps     = pm_ARTTV.GradDescSteps;
    recon             = pm_ARTTV.initial(:);
    alpha             = pm_ARTTV.alpha;
    epsilon           = pm_ARTTV.epsilon;
    sample_rate       = pm_ARTTV.sample_rate;
    
    % variables for monitoring
    visualize         = pm_ARTTV.visualize;
    disp_prog         = pm_ARTTV.disp_prog;
        
    % norms             = [];
    norms = zeros(iterations,1);
    errors = zeros(iterations,1);
    samples        = [];
    submetrics = ProblemSetup.empty_sub;
    submetric_keys = fieldnames(submetrics); % will be used for updating metrics
    % should be: {'t1': [], 't2':[], 't3':[]; 't4':[], 't5':[], 't6':[],
    %             't7':[], 't8':[]}
    overallmetrics    = [];
    
    if visualize
        fig_rows = 1;
        fig_cols = 2;
        figure(); set(gcf, 'units', 'normalized', 'outerposition', [0.25 0.25 .5 .5]);
        tiledlayout(fig_rows,fig_cols);
        
    end
    if disp_prog
        wb = waitbar(0, 'Iterating...'); 
    end
    for iter = 1:iterations
        if disp_prog
            waitbar((iter-1)/iterations, wb, sprintf('ART-TV Progress: %i of %i Iterations', ...
                                                 iter, iterations));
        end
        prev_reconstruction = recon;

        % (B) data projection steps
        for m = 1:NData
            proj_idx = proj_order(m);
            % basis vector
            Ai = A(proj_idx,:)';
            AiNorm = norm(Ai);
            if AiNorm>0
                p = projections(proj_idx);
                recon = recon + Ai*(p-Ai'*recon)/(AiNorm^2);
            end
        end
        
        % (C)-positivity step; calculate change (to initialize distance for grad
        %           desc step
        recon(recon<0)=0;

        if visualize
            nexttile(1);
            imshow(reshape(recon,[N,N]),[], 'InitialMagnification', 'fit');
            title('ARTTV-Recon post ART step')
        end
        % (D)-tv gradient descent initialization
        grad_desc_size = norm(recon-prev_reconstruction,2);
        % original exit condition
        % curr_dist = norm(reconstruction-prev_reconstruction,2);
        % % this is exit conditionused in L1dL2, tolerance 1e-5
        curr_dist = norm(recon-prev_reconstruction,'fro')/sqrt(numel(recon));
        
        % exit condition, occurs before we do gradient descent
        if curr_dist<exit_error
            [sub, overall] = localized_metric_v2(reshape(img,[N,N]), reshape(recon,[N,N]), ...
                                                 ProblemSetup);
            % 8 total fields, submetrics
            for i = 1:length(submetric_keys)
                submetrics.(submetric_keys{i}) = [submetrics.(submetric_keys{i}); ...
                                                    sub.(submetric_keys{i})];
            end                            
            % submetrics = cat(3, submetrics, sub);
            overallmetrics = [overallmetrics; overall];
            break
        end
        % (D')grad descent
        for m = 1:GradDescSteps
            v = dTV(recon,N,epsilon);
            v = v/norm(v);
            recon = recon-alpha*grad_desc_size*v;
        end
        norms(iter) = norm(recon-prev_reconstruction,2);
        errors(iter) = norm(recon - img,2);
        if visualize
            nexttile(1);
            imshow(reshape(recon,[N,N]),[], 'InitialMagnification', 'fit');
            title('ARTTV-Recon post TV step')
            nexttile(2);
            
            yyaxis left
            plot(norms(norms>0), 'b'); 
            ylabel('Exit Criteria (Change)');
        
            yyaxis right
            plot(errors(errors>0), 'r--'); 
            ylabel('Error');
        
            xlabel('Iterations');
            title('Exit Criteria and Error of Reconstruction by Iteration'); 
            pause(.0000001);
            
        end
        
        [sub, overall] = localized_metric_v2(reshape(img,[N,N]), reshape(recon,[N,N]), ...
                                             ProblemSetup);
        % 8 total fields, submetrics
        for i = 1:length(submetric_keys)
            submetrics.(submetric_keys{i}) = [submetrics.(submetric_keys{i}); sub.(submetric_keys{i})];
        end
        
        overallmetrics = [overallmetrics; overall];
        if mod(iter, sample_rate)==0
            samples = [samples; recon(:)'];
        end
    end
    metrics.sub = submetrics;
    metrics.full = overallmetrics;
    if disp_prog
        close(wb);
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