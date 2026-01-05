from argparse import ArgumentParser
import os

### Get command line arguments ###
def GetArgs():
    """Get command line arguments"""
    parser = ArgumentParser()

    parser.add_argument("-dir", "--directory", type=str, help="searches in subdirectories for ORCA output files")

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    args = GetArgs()

### note I only want this to contain key word arguments to dictate decisions - no other code







### uses command line arguments ###
base_dir = args.directory if args.directory else '.'

### finds directories ###
exclude = ('_SPC', 'IRC')
ts_dirs = [
    entry.name
    for entry in os.scandir(base_dir)
    if entry.is_dir()
    and entry.name.startswith('TS')
    and not any(x in entry.name for x in exclude)
]

print(ts_dirs)
### extract the data ###
for name in ts_dirs:
    sm_SPC = os.path.join(base_dir, name + '_IRC_B_SPC', name + '_IRC_B_SPC.out')
    p_SPC = os.path.join(base_dir, name + '_IRC_F_SPC', name + '_IRC_F_SPC.out')
    ts_SPC = os.path.join(base_dir, name + '_SPC', name + '_SPC.out')
    sm = os.path.join(base_dir, name + '_IRC_B', name + '_IRC_B.out')
    p = os.path.join(base_dir, name + '_IRC_F', name + '_IRC_F.out')
    ts = os.path.join(base_dir, name, name + '.out')
    ts_XYZ = os.path.join(base_dir, name, name + '.xyz')
    print(ts_XYZ)

