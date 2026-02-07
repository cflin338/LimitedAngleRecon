function hull = FindHull(ProblemSetup)
    hull = zeros(ProblemSetup.N, ProblemSetup.N);
    [x,y] = meshgrid(1:ProblemSetup.N, 1:ProblemSetup.N);
    cx = (ProblemSetup.N+1)/2;   % center
    cy = (ProblemSetup.N+1)/2;
    hull = ((x - cx).^2 + (y - cy).^2) <= (ProblemSetup.N/2)^2;
    hull = hull>0;
end