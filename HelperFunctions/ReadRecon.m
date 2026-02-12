function [recon,recon_slice] = ReadRecon(file_name, slice, XYdim)
    fid = fopen(file_name, 'r');
    recon = [];
    nextline = fgetl(fid);
    counter = 1;
    while ischar(nextline)
        row = sscanf(nextline, '%f')';  
        recon = [recon; row]; 
        nextline = fgetl(fid); 
        counter = counter + 1;
    end
    recon_slice = recon((slice-1)*XYdim+1:slice*XYdim,:);
end