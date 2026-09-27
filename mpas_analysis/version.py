import re

# must match the git tag exactly, including any rc suffix
__version__ = '1.16.0rc1'

# the numeric part of the version, e.g. (1, 16, 0) for '1.16.0rc1'
__version_info__ = tuple(
    int(part)
    for part in re.split(r'[^\d.]', __version__, maxsplit=1)[0].split('.')
)
