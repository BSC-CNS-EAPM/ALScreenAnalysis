"""Shared validation and manifest loading for analogue analyses."""
from pathlib import Path
import pandas as pd


def read_table(path, required, numeric=(), unique=None):
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f'{path}: missing columns {sorted(missing)}')
    for col in required:
        if df[col].str.strip().eq('').any():
            raise ValueError(f'{path}: blank values in {col}')
    for col in numeric:
        df[col] = pd.to_numeric(df[col], errors='raise')
        import numpy as np
        if not np.isfinite(df[col]).all():
            raise ValueError(f'{path}: nonfinite values in {col}')
    if unique and df.duplicated(unique).any():
        raise ValueError(f'{path}: duplicate keys {unique}; select best poses first')
    return df


def manifest(path, extra=()):
    path = Path(path).resolve()
    required = ['target', 'hit', 'constraint', 'analogue_csv', 'hit_csv', 'pele_csv', *extra]
    df = read_table(path, required, unique=['target', 'hit', 'constraint'])
    if df.empty:
        raise ValueError('Manifest must contain at least one row')
    for col in ['analogue_csv', 'hit_csv', 'pele_csv', 'visual_dir', 'pose_dir']:
        if col in df:
            df[col] = df[col].map(lambda v: str((path.parent / v).resolve()) if v else '')
    # These fields form output paths. IDs inside data files need no filename parsing.
    for col in ['target', 'hit', 'constraint']:
        for value in df[col]:
            if value in ('.', '..') or '/' in value or '\\' in value:
                raise ValueError(f'Unsafe path component in {col}: {value!r}')
    return df


def suffix(row, kind='pele'):
    return row.get(f'{kind}_suffix', '') or row['constraint']


def molecule_name(identifier, label):
    return f'{identifier}_{label}'


def single_hit(df, hit, id_col, context):
    selected = df.loc[df[id_col].eq(hit)]
    if len(selected) != 1:
        raise ValueError(f'{context}: expected one row for hit {hit!r}, found {len(selected)}')
    return selected.iloc[0]


def optional_number(row, key):
    value = row.get(key, '')
    if value == '':
        return None
    import math
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f'{key} must be finite')
    return result


def output_stem(value):
    path = Path(value)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path
