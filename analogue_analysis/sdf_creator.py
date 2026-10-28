"""Batch CSV-to-SDF conversion through an external helper."""
from batch_common import parser, files, run


def main():
    p = parser(__doc__)
    p.add_argument('--id-column', default='name')
    p.add_argument('--smiles-column', default='smiles')
    p.add_argument('--properties', nargs='+', default=['space', 'smiles'])
    p.add_argument('--delimiter', default=',')
    a = p.parse_args()
    for directory in a.directories:
        for path in files(directory, '*.csv'):
            run(a, ['--csv', path, '--id', a.id_column, '--smi', a.smiles_column,
                    '--params', *a.properties, '-d', a.delimiter, '-o', path.with_suffix('')])


if __name__ == '__main__':
    main()
