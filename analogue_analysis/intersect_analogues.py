"""Intersect matching MolDB filenames in two researcher-selected folder trees."""
from pathlib import Path
from batch_common import parser, files, run


def main():
    p = parser(__doc__)
    p.add_argument('--other-base', required=True, help='Contains counterpart folders named like each input folder')
    p.add_argument('--output-dir', required=True)
    p.add_argument('--similarity', type=float, default=1.0)
    p.add_argument('--no-sdf', action='store_true')
    a = p.parse_args()
    if not 0 <= a.similarity <= 1:
        p.error('Similarity must lie between 0 and 1')
    seen = set()
    for directory in a.directories:
        name = Path(directory).resolve().name
        if name in seen:
            p.error('Input folders must have distinct basenames to avoid output collisions')
        seen.add(name)
        counterpart = Path(a.other_base).resolve() / name
        if not counterpart.is_dir():
            raise NotADirectoryError(counterpart)
        matched = 0
        for path in files(directory, '*.moldb'):
            other = counterpart / path.name
            if not other.is_file():
                print(f'Skipping unmatched file: {path}')
                continue
            matched += 1
            outdir = Path(a.output_dir).resolve() / name
            if not a.dry_run:
                outdir.mkdir(parents=True, exist_ok=True)
            run(a, ['--moldb1', path, '--moldb2', other,
                    '-o', outdir / (path.stem + '_intersect'), '--sim', a.similarity,
                    *([] if a.no_sdf else ['--sdf'])])
        if not matched:
            raise ValueError(f'No matching MolDB filenames for {directory}')


if __name__ == '__main__':
    main()
