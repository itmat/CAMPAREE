#!/usr/bin/env python3

import sys

sp_dir = sys.argv[2].rstrip()
id     = sys.argv[3]
in_dir = sys.argv[4]
end1   = sys.argv[5]
end2   = sys.argv[6]
resources = 0
output    = 0
with open(sys.argv[1]) as in_yaml, open(sys.argv[7], 'w') as out_yaml:
    for line in in_yaml:
        line = line.rstrip();
        if line == 'resources:':
            resources = 1
            out_yaml.write(line + '\n')
        elif line == '    directory_path:':
            if resources:
                out_yaml.write(line + ' ' + sp_dir + '/resources/\n')
            elif output:
                out_yaml.write(line + ' ' + sp_dir + '/results/\n')
        elif line == '    fastq_directory_path:':
            out_yaml.write(line + ' ' + sp_dir + '/' + in_dir + '/\n')
        elif line == '    data:':
            out_yaml.write(line + '\n')
            out_yaml.write('        \'' + id + '\':\n')
            out_yaml.write('            fastq_files:\n')
            out_yaml.write('                - ' + end1 + '\n')
            out_yaml.write('                - ' + end2 + '\n')
            out_yaml.write('            pooled: False\n')
            out_yaml.write('            gender: Male\n')
            out_yaml.write('            molecule_count: 1000000\n')
            out_yaml.write('            optional_inputs:\n')
        elif line == 'output:':
            resources = 0
            output    = 1
            out_yaml.write(line + '\n')
        else:
            out_yaml.write(line + '\n')