"""Write a simple F90 module that defines PF_VERSION and PF_GIT_VERSION."""

import os
import re
import sys


# config

version = 'v0.1.0'
root    = os.environ.get('PFASST', '.')


# checks

if len(sys.argv) != 2:
    print( "ERROR: must pass output file name." )
    raise SystemExit

# get the git version

try:
    git_head_file = os.path.join(root, '.git', 'HEAD')
    f = open(git_head_file)
    m = re.match(r'ref: (.+)', f.readline())
    ref = m.group(1)
    f.close()

    git_head_file = os.path.join(root, '.git', ref)
    f = open(git_head_file)
    git_version = f.readline().rstrip()
    f.close()

except:
    git_version = 'not_available'


# define f90 module

src = '''module pf_mod_version
  character(len=*), parameter :: PF_VERSION = '{version}'
  character(len=*), parameter :: PF_GIT_VERSION = '{git_version}'
end module pf_mod_version
'''

src = src.format(version=version, git_version=git_version)

# if the output file already exists, check to see if anything has changed

output_file = sys.argv[1]

if os.path.exists(output_file):
    with open(output_file, 'r') as f:
        existing_src = f.read()

    if src == existing_src:
        raise SystemExit

# write...

with open(output_file, 'w') as f:
    f.write(src)
