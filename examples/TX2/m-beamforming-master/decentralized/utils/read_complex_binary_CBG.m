function v = read_complex_binary_CBG(fid, count)

    position_before = ftell(fid);
    t = fread(fid, [2, count], 'float');
    position_after = ftell(fid);
    lengthMyRead = length(t);
    
    if lengthMyRead==count
        v = t(1,:) + t(2,:)*1i;
        [r, c] = size (v);
        v = reshape (v, c, r);
    elseif lengthMyRead>0
        fseek(fid,position_before - position_after,'cof');
%         fprintf('Correcting cursor for misalignment by: %d\n',position_before - position_after);
        fprintf('LOG - Correcting file Misalignment\n');
        v = [];
    else
        v = [];
    end
    
    
% EOF
