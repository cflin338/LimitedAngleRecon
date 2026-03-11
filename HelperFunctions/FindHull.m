function hull = FindHull(ProblemSetup, rad)
    hull = zeros(ProblemSetup.N, ProblemSetup.N);
    [x,y] = meshgrid(1:ProblemSetup.N, 1:ProblemSetup.N);
    cx = (ProblemSetup.N+1)/2;   % center
    cy = (ProblemSetup.N+1)/2;
    hull = ((x - cx).^2 + (y - cy).^2) <= rad^2;
    hull = hull>0;
end