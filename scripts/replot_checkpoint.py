"""Re-render stage 2/3 figures without rerunning clustering or changing h5ad.

Example: python replot_checkpoint.py --config config.json --checkpoint final.h5ad
         --groupby subtype --output figure_revision
Output must be a fresh/empty directory, keeping previous figures intact.
"""
import argparse
import json
import os
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--checkpoint', required=True)
    parser.add_argument('--groupby', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    destination = Path(args.output)
    if destination.exists() and any(destination.iterdir()):
        raise ValueError('Use a fresh output directory for figure-only revisions')
    cache = destination.parent / '.plot_runtime_cache'
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault('NUMBA_CACHE_DIR', str(cache / 'numba'))
    os.environ.setdefault('MPLCONFIGDIR', str(cache / 'matplotlib'))
    dll = Path(sys.prefix) / 'Library' / 'bin'
    if dll.exists():
        os.environ['PATH'] = str(dll) + os.pathsep + os.environ.get('PATH', '')
        if hasattr(os, 'add_dll_directory'):
            dll_handle = os.add_dll_directory(str(dll))
    import scanpy as sc
    from human_review import draw_review
    config = json.loads(Path(args.config).read_text(encoding='utf-8-sig'))
    adata = sc.read_h5ad(args.checkpoint)
    if args.groupby not in adata.obs:
        raise ValueError(f'Unknown groupby: {args.groupby}')
    adata.obs[args.groupby] = adata.obs[args.groupby].astype('category')
    draw_review(adata, args.groupby, destination, config)
    sys.stdout.write(f'Rendered {len(list(destination.glob("*.png")))} PNG/PDF pairs to {destination}\n')


if __name__ == '__main__':
    main()
