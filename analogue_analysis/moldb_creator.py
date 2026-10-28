"""Batch SDF-to-MolDB conversion through an external helper."""
from batch_common import parser, files, run


def main():
    a = parser(__doc__).parse_args()
    for directory in a.directories:
        inputs = files(directory, '*.sdf')
        run(a, ['--sdfs', *inputs, '-o', *(f.with_suffix('') for f in inputs)])


if __name__ == '__main__':
    main()
