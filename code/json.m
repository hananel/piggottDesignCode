function ret = json(filename, s);

% Convert s to a json formatted string, writing it to filename
% (We could make a version that wrote to a string if required)
% s is a struct with fields of strings (can you have something else?)
% So it is converted to json like:
%
% s = struct()
% s.a = 1.2
% s.b = 3
% s.c = 'asd'
% s.d = struct() # not supported yet
% s.e = [1,2,3]  # not supported yet
% 
% json('output.txt', s)
% $ cat output.txt
%
% {
%  "a": 1.2,
%  "b": 3,
%  "c": "asd"
% }\n
%

fd = fopen(filename, 'w+');
fprintf(fd, '{\n');
fields = fieldnames(s);
n = size(fields)(1);
for i = 1:n
    name = fields{i};
    value = getfield(s, name);
    fprintf(fd, ' "%s": ', name);
    if (ischar(value))
        fprintf(fd, '"%s"', value)
    elseif (int64(value) == value)
        fprintf(fd, '%d', value);
    elseif (isfloat(value))
        fprintf(fd, '%f', value);
    end
    fprintf(fd, ',\n');
end
fprintf(fd, '}\n');
fclose(fd)
