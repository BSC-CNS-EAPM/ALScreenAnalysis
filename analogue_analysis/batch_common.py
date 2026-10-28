"""Shared safe invocation for external compound-analysis helper scripts."""
import argparse
from pathlib import Path
import shlex
import subprocess
import sys


def parser(description):
    p = argparse.ArgumentParser(description=description)
    p.add_argument('--directories', nargs='+', required=True, help='Input folders; no implicit target list')
    p.add_argument('--script', required=True, help='External helper script path')
    p.add_argument('--python', default=sys.executable, help='Interpreter for external helper')
    p.add_argument('--dry-run', action='store_true', help='Print commands without executing')
    return p


def run(args, extra):
    helper = Path(args.script).resolve()
    if not helper.is_file():
        raise FileNotFoundError(helper)
    command = [args.python, str(helper), *map(str, extra)]
    print(shlex.join(command), flush=True)
    if not args.dry_run:
        subprocess.run(command, check=True)


def files(directory, pattern):
    path = Path(directory).resolve()
    if not path.is_dir():
        raise NotADirectoryError(path)
    result = sorted(f for f in path.glob(pattern) if f.is_file())
    if not result:
        raise ValueError(f'No {pattern} files in {path}')
    return result
