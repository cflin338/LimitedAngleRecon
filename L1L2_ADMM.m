function [u,metrics, samples] = L1L2_ADMM(ProblemSetup, pm_L1dL2)
    % Limited-angle CT reconstruction via the L1/L2 minimization∗
    %       Rahimi 2019
    % Limited-Angle CT Reconstruction via the $L_1/L_2$ Minimization
    %       Wang 2021
    %
    % implemented unconstrained, bounded box version
    % L1/L2 on gradient in unconstrained formulation
    % -----min u st norm(grad(u),1)/norm(grad(u),2) + lambda/2*norm(A*u-f,2)^2
    % ADMM:
    % -----min u st norm(grad(u),1)/norm(h,2) + lambda/2*norm(A*u-f,2)^2 st 
    %                                                               h = grad(u)                       
    % augmented lagrangian:
    % -----L = norm(grad(u),1)/norm(h,2) +  lambda/2*norm(A*u-f,2)^2 + 
    %                   dot(rho2*b2, grad(u)-h) + rho2/2*norm(h-grad(u),2)^2
    % update: solve u, solve h
    %           h = tau*g or e
    %           b2 = b2 + grad(u) - h
    % updating u, admm u sub problem, introduce d var
    
    % General Variables
    A                 = ProblemSetup.A;
    projections       = ProblemSetup.projections;
    img               = ProblemSetup.img;
    N                 = ProblemSetup.N;
    sample_rate       = pm_L1dL2.sample_rate;
    % sample_rate       = 1;
    % variables for monitoring
    visualize         = pm_L1dL2.visualize;
    disp_prog         = pm_L1dL2.disp_prog;
    
    % Method Specific Variables
    rho1              = pm_L1dL2.rho1; 
    rho2              = pm_L1dL2.rho2; 
    gamma             = rho1+rho2;
    beta              = pm_L1dL2.beta; 
    lambda            = pm_L1dL2.lambda; 
    iterations        = pm_L1dL2.maxit; 
    subiterations     = pm_L1dL2.Imaxit; 
    inner_epsilon     = pm_L1dL2.Itol; 
    outer_epsilon     = pm_L1dL2.tol;
    upper_bound       = pm_L1dL2.upper_bound;
    lower_bound       = pm_L1dL2.lower_bound;
    
    samples        = [];

    % initial guess
    u = reshape(pm_L1dL2.initial,[N,N]);
    e = zeros([N,N]); 
    
    dx = zeros([N,N]);
    dy = zeros([N,N]);
    bx = zeros([N,N]);
    by = zeros([N,N]);
    
    hx = zeros([N,N]);
    hy = zeros([N,N]);
    cx = zeros([N,N]);
    cy = zeros([N,N]);
    
    if visualize
        fig_rows = 1;
        fig_cols = 2;
        figure(5); set(gcf, 'units', 'normalized', 'outerposition', [0.25 0.25 .5 .5]);
        tiledlayout(fig_rows,fig_cols);
    end
    
    errors            = [];

    submetrics = ProblemSetup.empty_sub;
    submetric_keys = fieldnames(submetrics); % will be used for updating metrics
    overallmetrics    = [];
    
    for iter = 1:iterations
        if disp_prog
            fprintf("%i of %i iterations \n",iter,iterations); 
        end
        prev_u = u;
        % update u
        % introduce v for subproblem
        v = u;
        for j = 1:subiterations
            prev_subu = u;
            rhs = (lambda * reshape(A' * projections,[N,N]) + ...
                    rho1 * Dxt(dx-bx) + rho1 * Dyt(dy-by) + ...
                    rho2 * Dxt(hx-cx) + rho2 * Dyt(hy-cy) + ...
                    beta * (v - e));

            [u,~] = conjgrad_b(A,rhs(:),u(:),N,N,gamma,lambda,beta);
            u = reshape(u, [N,N]);
            
            % update gradient of u; Dxu/Dyu
            Dxu = Dx(u);
            Dyu = Dy(u);

            % update d; update dx, dy
            hnorm = sqrt(norm(hx(:))^2+norm(hy(:))^2); % norm of h overall
            dx = shrink(Dxu+bx, 1/(rho1*hnorm));
            dy = shrink(Dyu+by, 1/(rho1*hnorm));
            % update b; update bx, by; this is the c in previous code/ the
            %                           paper
            bx = bx + (Dxu - dx);
            by = by + (Dyu - dy);

            % update v
            % bounding within [0,1] if we assume phantom used
            v = u+e;
            v(v<lower_bound)=lower_bound; 
            v(v>upper_bound)=upper_bound;
            % update e in subloop if box constrained; no change
            e = e + u - v;

            subu_stepsize = norm(prev_subu-u,'fro') / norm(u,'fro');
            
            if subu_stepsize<inner_epsilon
                break
            end
        end

        % update h, using updated u and updated h
        d1 = Dxu + cx;
        d2 = Dyu + cy;
        etha = sqrt(norm(d1(:))^2+norm(d2(:))^2); 
        c = norm(dx(:),1)+norm(dy(:),1);
        % update e norm size, only need to do this if performing unboxed
        %   iteration
        if etha == 0
            hx = (c/rho)^(1/3)*ones(size(d1))/sqrt(numel(d1)*2);
            hy = (c/rho)^(1/3)*ones(size(d1))/sqrt(numel(d1)*2);
        else
            a = 27*c/(rho2*(etha^3)) + 2;
            C = ((a + (a^2 - 4)^0.5)/2)^(1/3);
            tau = (1 + C + 1/C)/3;
            hx = tau*d1;
            hy = tau*d2;
        end
        % update c-cx/cy, b2 in previous code
        cx = cx + Dxu - hx;
        cy = cy + Dyu - hy;

        stepsize = norm(u-prev_u,'fro')/sqrt(numel(u));
        errors = [errors, stepsize];
        [sub, overall] = localized_metric_v2(reshape(img,[N,N]), reshape(u,[N,N]), ...
                                             ProblemSetup);
        % 8 total fields, submetrics
        for metric_idx = 1:length(submetric_keys)
            submetrics.(submetric_keys{metric_idx}) = [submetrics.(submetric_keys{metric_idx}); ...
                                                sub.(submetric_keys{metric_idx})];
        end  

        overallmetrics = [overallmetrics; overall];
        
        if visualize
            nexttile(1);
            imshow(u,[], 'InitialMagnification', 'fit');
            title('L1/L2 Reconstruction');
            nexttile(2);
            plot(errors);
            title('Exit Criteria'); xlabel('Iterations'); ylabel('Frobenius Norm');
            pause(.00001);
        end
        
        if stepsize < outer_epsilon
            break
        end
        if mod(iter, sample_rate)==0
            samples = [samples; u(:)'];
        end
    end
    metrics.sub = submetrics;
    metrics.full = overallmetrics;
    
end

function v = shrink(s,lambda)
    v = sign(s).*(max(abs(s) - lambda , 0));
end

function d = Dx(u)
    ux = size(u,1);
    d = zeros(ux,ux);
    d(:,2:ux) = u(:,2:ux)-u(:,1:ux-1);
    d(:,1) = u(:,1)-u(:,ux);
end

function d = Dxt(u)
    ux = size(u,1);
    d = zeros(ux,ux);
    d(:,1:ux-1) = u(:,1:ux-1)-u(:,2:ux);
    d(:,ux) = u(:,ux)-u(:,1);
end

function d = Dy(u)
    ux = size(u,1);
    d = zeros(ux,ux);
    d(2:ux,:) = u(2:ux,:)-u(1:ux-1,:);
    d(1,:) = u(1,:)-u(ux,:);

end

function d = Dyt(u)
    ux = size(u,1);
    d = zeros(ux,ux);
    d(1:ux-1,:) = u(1:ux-1,:)-u(2:ux,:);
    d(ux,:) = u(ux,:)-u(1,:);

end

function [u,rsnew] = conjgrad_b(R,rhs,u,rows,cols,gamma,lambda,beta)
% https://github.com/wangcmath/limited_angle_CT_L1dL2/blob/main/mCTrecon_L1dL2_unconst.m
% output x is the entire update: rhs / (...)
%               lambda * AtA - gamma*laplacian+ beta*I
% 
    % (A=A, b=rhs(:), x=u(:), pm=cgpm)
    r=rhs-Ax_b(R,u,rows,cols,gamma,lambda,beta);
    % 
    p=r;
    rsold=r'*r;
 
    for i=1:sqrt(rows)
        Ap=Ax_b(R,p,rows,cols,gamma,lambda,beta);
        alpha=rsold/(p'*Ap);
        u=u+alpha*p;
        r=r-alpha*Ap;
        rsnew=r'*r;
        if sqrt(rsnew)<1e-2
              break;
        end
        p=r+rsnew/rsold*p;
        rsold=rsnew;
    end
end
    
function y = Ax_b(R,u, rows,cols,gamma,lambda,beta)
% (P=A, x=u(:), pm=cgpm)
g = fspecial('laplacian',0);
y = imfilter(reshape(u,rows,cols),g, 'circular');
y =  lambda*R'*(R*u(:))-gamma*y(:)+beta*u(:);
end