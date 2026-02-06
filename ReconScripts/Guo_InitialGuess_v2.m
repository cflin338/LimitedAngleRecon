% Iterative Image Reconstruction for Limited-Angle CT Using Optimized Initial Image
%      Jingyu Guo 2016

function initial = Guo_InitialGuess_v2(ProblemSetup,J,K, visualize)
    A                 = ProblemSetup.A;
    projections       = ProblemSetup.projections;
    img               = ProblemSetup.img;
    bins              = ProblemSetup.bins;
    angles            = ProblemSetup.angles;
    N                 = ProblemSetup.N;
    projections = reshape(projections, [bins,length(angles)]);

    if nargin==3
        visualize=false;
    end

    initial = iradon(projections, angles);
    initial = imresize(initial,N/size(initial,1));

    if visualize
        figure; tiledlayout(2,2);
        nexttile(1); imshow(img,[]);
        nexttile(2); imshow(initial,[]);
    end
    
    L_contour = [];
    R_contour = [];
        
    base = sparse(N,N);
    projections = projections(:);
    for p = 1:length(projections)
        if projections(p)==0
            base = base | reshape(A(p,:),[N,N]);
        end
    end
    
    initial = initial.*(1-base);
    if visualize
        nexttile(3);
        imshow(full(base));
    end       
    for row = 1:size(initial,1)
        % find contours, append to L, R
        % this is first and last nonzero of projections
        delta = find(diff(base(row,:)));
        if isempty(delta)
            L_contour = [L_contour, 0];
            R_contour = [R_contour, 0];
        else
            left_nonzero = delta(1)+1;
            L_contour = [L_contour, left_nonzero];
            right_nonzero = delta(2);
            R_contour = [R_contour, right_nonzero];
        end
    end
    % find top of contour
    found = false;
    m_top = 1;
    while ~found
        delta = find(diff(base(m_top,:)));
        if isempty(delta)
            m_top=m_top+1;
        else
            found=true;
        end
    end
    % find bot of contour
    found = false;
    m_bot = N;
    while ~found
        delta = find(diff(base(m_bot,:)));
        if isempty(delta)
            m_bot=m_bot-1;
        else
            found=true;
        end
    end
    
    % calculate center of symmetry for top J+1 rows
    S_top = 0;
    for j = 0:J
        S_top = S_top + (R_contour(m_top+j)+L_contour(m_top+j))/2;
    end
    S_top = round(S_top / (J+1));
    
    S_bot = 0;
    for j = 0:J
        S_bot = S_bot + (R_contour(m_bot-j)+L_contour(m_bot-j))/2;
    end
    S_bot = round(S_bot / (J+1));
    
    Centerpoint = round((S_top+S_bot)/2);
    S_top = Centerpoint;
    S_bot = Centerpoint;
    % loop from m+J to N/2
    for row = m_top+J:N/2
        % reflect right contour to left, around S_top center
        L_contour(row) = max(1,S_top - (R_contour(row)-S_top));
        initial(row, L_contour(row):L_contour(row)+min(K,S_top-L_contour(row))) = ...
            flip(initial(row, R_contour(row)-min(K,S_top-L_contour(row)):R_contour(row)));
        initial(row,1:L_contour(row)-1)=0;
        initial(row,R_contour(row)+1:N)=0;
    end
    
    % loop from m+J to N/2
    for row = m_bot-J:-1:N/2+1
        R_contour(row) = min(S_bot + (S_bot-L_contour(row)),N);
        initial(row, R_contour(row)-min(K,R_contour(row)-S_bot):R_contour(row)) = ...
            flip(initial(row, L_contour(row):L_contour(row)+min(K,R_contour(row)-S_bot)));
        initial(row,R_contour(row)+1:N)=0;
        initial(row,1:L_contour(row)-1)=0;
    end

    if visualize
        nexttile(4);
        imshow(initial);
    end
    initial=initial(:);
end

