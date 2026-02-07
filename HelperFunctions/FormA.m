function A = FormA(ProblemSetup)

    det_width = ProblemSetup.ssd_tsize;
    det_bins = ProblemSetup.tbins;
    ang_count = ProblemSetup.angle_count;
    N = ProblemSetup.N;
    
    rows = [];
    cols = [];
    vals = [];

    dt = ProblemSetup.tbinsize;
    angles = linspace(0,2*pi,ProblemSetup.angle_count+1);
        angles(end)=[];
    R = ProblemSetup.source_radius; %assume distance from source to detector is 2R
    dx = ProblemSetup.voxel_width;
    
    % assume center is 0,0
    % assume starting positions: [0, -R], [0, R]
    detector_base_locations = [linspace(-det_width/2+dt/2,det_width/2-dt/2,det_bins)', R*ones(det_bins,1)];
    orig_src_loc = [0, -R];
    vx = linspace(-N/2*dx+dx, N/2*dx-dx, N);
    vy = linspace(-N/2*dx+dx, N/2*dx-dx, N);
    % ang_idx = 1; bin_angle = angles(ang_idx);
    
    for ang_idx = 1:ang_count
        bin_angle = angles(ang_idx);
        sinA = sin(bin_angle);
        cosA = cos(bin_angle);

        for det_idx = 1:det_bins
            % cum_val = sinogram(ang_idx, det_idx);
            row_id = det_idx + (ang_idx-1)*det_bins;
            % determine start, stop location
            det_loc_orig = detector_base_locations(det_idx,:);
            det_loc = [det_loc_orig(1)*cosA - det_loc_orig(2)*sinA, ...
                            det_loc_orig(1)*sinA + det_loc_orig(2)*cosA];
            source_loc = [orig_src_loc(1)*cosA - orig_src_loc(2)*sinA, ...
                            orig_src_loc(1)*sinA + orig_src_loc(2)*cosA];
            % use siddon method to calculate the pixels
            [pix_idx, weights] = siddon2D(det_loc, source_loc, vx, vy);
            rows = [rows; repmat(row_id,length(pix_idx),1)];
            cols = [cols; pix_idx(:)];
            vals = [vals; weights(:)];

        end
    end
    A = sparse(rows, cols, vals, det_bins*ang_count, N*N);
    
end


function [indices,lengths] = siddon2D(det_loc, source_loc, xv, yv)
sx = det_loc(1); sy = det_loc(2);
dx= source_loc(1); dy = source_loc(2);
% Exact 2D Siddon ray–pixel intersection
% xv, yv: voxel centers
% Returns:
%   indices -> linear indices
%   lengths -> intersection lengths in pixels

dxv = xv(2)-xv(1);
dyv = yv(2)-yv(1);
N = length(xv);

xmin = xv(1)-dxv/2; xmax = xv(end)+dxv/2;
ymin = yv(1)-dyv/2; ymax = yv(end)+dyv/2;

dir = [dx-sx, dy-sy];
L = norm(dir);
dir = dir/L;

% intersection t for box boundaries
tx = ([xmin xmax]-sx)/dir(1);
ty = ([ymin ymax]-sy)/dir(2);

t0 = max(min(tx), min(ty));
t1 = min(max(tx), max(ty));

if t1 <= t0
    indices = [];
    lengths = [];
    return
end

% start point inside voxel
p = [sx,sy] + t0*dir;
ix = floor((p(1)-xmin)/dxv)+1;
iy = floor((p(2)-ymin)/dyv)+1;

ix = min(max(ix,1),N);
iy = min(max(iy,1),N);

% step direction
stepx = sign(dir(1)); stepy = sign(dir(2));
if stepx==0, txnext=inf; else txnext = t0 + dxv/abs(dir(1)); end
if stepy==0, tynext=inf; else tynext = t0 + dyv/abs(dir(2)); end

t = t0;
indices = [];
lengths = [];

while t < t1 && ix>=1 && ix<=N && iy>=1 && iy<=N
    % find next crossing
    if txnext <= tynext
        tnew = txnext;
        axis = 1;
    else
        tnew = tynext;
        axis = 2;
    end

    seglen = min(tnew,t1) - t;

    if seglen > 0
        indices(end+1,1) = sub2ind([N,N], ix, iy);
        lengths(end+1,1) = seglen*L;
    end

    % step index
    if axis==1
        ix = ix + stepx;
        txnext = txnext + dxv/abs(dir(1));
    else
        iy = iy + stepy;
        tynext = tynext + dyv/abs(dir(2));
    end

    t = tnew;
end

end
