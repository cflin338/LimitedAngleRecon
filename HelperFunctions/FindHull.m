function hull = FindHull(sinogram, A, threshold)
    sinogram = sinogram == 0;
    % sum up all projection vectors that miss the object
    miss_counts = sinogram*A;
    hull = miss_counts>threshold;
end