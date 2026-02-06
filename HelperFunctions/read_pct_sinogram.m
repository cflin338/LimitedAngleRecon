function [sino_raw, sino_slice] = read_pct_sinogram(file_name, vsize, tsize, angle_count, slice)
    fid = fopen(file_name, 'r');
    sino_raw = [];
    nextline = fgetl(fid);
    counter = 1;
    while ischar(nextline)
        row = sscanf(nextline, '%f')';  
        sino_raw = [sino_raw; row]; 
        nextline = fgetl(fid); 
        counter = counter + 1;
    end

    fclose(fid);
    sino_slice = sino_raw((slice-1)*angle_count+1:slice*angle_count, :);
    % reformat so each groupings instead correspond with angles
    % reformatted_sinogram = zeros(angle_count, vsize, tsize);
    % for angle = 1:angle_count
    %     detector = zeros(vsize, tsize);
    %     for v = 1:vsize
    %         detector(v,:) = sino_raw( (v-1)*angle_count + angle ,:);
    %     end
    %     reformatted_sinogram(angle,:,:) = detector;
    % end
end