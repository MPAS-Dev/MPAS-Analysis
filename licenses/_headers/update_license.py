#!/usr/bin/env python
import fileinput
import glob

extensions = ['py', 'sh', 'bash', 'html', 'css', 'js', 'ncl']

replace = {'Copyright (c) 2020 Triad National Security, LLC.':
           'Copyright (c) 2020 Triad National Security, LLC.',
           'Copyright (c) 2020 Lawrence Livermore National Security, LLC.':
           'Copyright (c) 2020 Lawrence Livermore National Security, LLC.',
           'Copyright (c) 2020 UT-Battelle, LLC.':
           'Copyright (c) 2020 UT-Battelle, LLC.'}

files = []
for ext in extensions:
    files.extend(glob.glob('./**/*.{}'.format(ext), recursive=True))

files.extend(glob.glob('./**/header-*', recursive=True))

for filename in files:
    print(filename)
    with fileinput.FileInput(filename, inplace=True, backup=None) as file:
        for line in file:
            for old in replace:
                new = replace[old]
                line = line.replace(old, new)
            print(line, end='')
