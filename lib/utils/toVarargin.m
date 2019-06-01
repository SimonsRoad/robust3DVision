function output = toVarargin(data)

argnames = fieldnames(data)';
output = cell(1, length(argnames)*2);
i = 1;
for args = argnames
   output{i} = args{1};
   output{i+1} = data.(args{1});
   i = i+2;
end

end