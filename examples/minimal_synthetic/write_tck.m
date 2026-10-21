function write_tck(fname, streamlines)
% WRITE_TCK  Write streamlines (cell array of 3 x N world-mm vertices) as an MRtrix .tck file.
hdr = sprintf(['mrtrix tracks\ndatatype: Float32LE\ncount: %d\n' ...
               'file: . %%d\nEND\n\n'], numel(streamlines));
off = numel(sprintf(hdr, 0)) + 8;                 % room for the offset digits
off = numel(sprintf(hdr, off));
txt = sprintf(hdr, off);
fid = fopen(fname, 'w', 'ieee-le');
fwrite(fid, txt, 'char');
fwrite(fid, zeros(1, off - numel(txt)), 'uint8');
for i = 1:numel(streamlines)
    fwrite(fid, single(streamlines{i}), 'float32');
    fwrite(fid, single(nan(3,1)), 'float32');
end
fwrite(fid, single(inf(3,1)), 'float32');
fclose(fid);
end
