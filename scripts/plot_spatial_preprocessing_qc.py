#!/usr/bin/env python3
"""Render source-linked donor QC only; never changes a spot mask or fits a model."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, StrMethodFormatter

DONORS = ['P1', 'P5', 'P8', 'P9', 'P10', 'P11', 'P12', 'P13']


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def flag(s):
    assert s.astype(str).str.lower().isin(['true', 'false']).all()
    return s.astype(str).str.lower().eq('true').to_numpy()


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument('--root', required=True, type=Path)
    a.add_argument('--output-dir', default='figures/spatial_preprocessing_QC')
    x = a.parse_args()
    root = x.root.resolve()
    assert root.name == 'IBS_A_OFFICIAL_RERUN_20261105'
    gate_path = root / 'manifest/PHASE1_GATE.json'
    assert json.loads(gate_path.read_text())['status'] == 'PASS'
    status_path = root / 'data/metadata/SPATIAL_PREPROCESS_ALL_DONORS_STATUS.json'
    assert json.loads(status_path.read_text())['status'] == 'PREPROCESS_COMPLETE'
    out = (root / x.output_dir).resolve()
    assert (root / 'figures').resolve() in out.parents, 'Output must remain in new figures directory'
    out.mkdir(parents=True, exist_ok=False)
    source_dir = out / 'source_data'
    source_dir.mkdir()
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'],
                         'font.size': 8, 'axes.spines.right': False, 'axes.spines.top': False,
                         'axes.linewidth': .7, 'legend.frameon': False, 'pdf.fonttype': 42,
                         'svg.fonttype': 'none'})
    contract = {
        'claim': 'Audit the actual importer/background/QC/manual-exclusion outcome separately for each donor; no disease association or model robustness is inferred.',
        'archetype': 'quantitative grid', 'backend': 'Python matplotlib',
        'biological_unit': 'one surgical reference donor per page; spots are descriptive observations, not independent donors',
        'statistics': 'Descriptive counts, empirical distributions and frozen thresholds only; no new P-values or significance tests.',
        'panels': ['all original UTT pixel coordinates and unchanged raw UMI', 'retained primary pixel coordinates and corrected expression',
                   'all original UTT source-defined z scores and fixed main/strict gates', 'original gene-row background distribution',
                   'raw/corrected detected-gene distributions', 'raw/corrected mitochondrial fractions'],
        'coordinate_integrity': 'Full-resolution pixel coordinates used only for these display maps; y-axis inverted. Not Visium array adjacency and no tissue-image scale bar claim.',
        'export': {'inches': [10.8, 7.3], 'PNG_dpi': 200, 'PDF_editable_text': True},
        'review_risks': ['Independent Python importer reimplementation, not MATLAB execution or full zonation.',
                         'Main and strict sensitivity masks are frozen upstream; this script never selects spots.',
                         'Strict-mask model sensitivity is not performed by these plots.', 'Raw integer and corrected fractional counts remain distinct.']
    }
    (out / 'FIGURE_CONTRACT.json').write_text(json.dumps(contract, indent=2) + '\n')
    sources = {str(p.relative_to(root)): sha(p) for p in [gate_path, status_path]}
    summary = []
    for donor in DONORS:
        folder = root / 'data/processed/spatial' / donor
        paths = [folder / 'UTT_QC_all_spots.tsv', folder / 'background_by_original_gene_row.tsv', folder / 'PREPROCESS_STATUS.json']
        sources.update({str(p.relative_to(root)): sha(p) for p in paths})
        qc = pd.read_csv(paths[0], sep='\t')
        bg = pd.read_csv(paths[1], sep='\t')
        st = json.loads(paths[2].read_text())
        keep = flag(qc.entered_primary_analysis)
        strict = flag(qc.entered_strict_z3_5_sensitivity)
        manual = flag(qc.in_author_manual_exclusion)
        assert int(keep.sum()) == st['n_retained'] and not (keep & manual).any()
        assert not (strict & ~keep).any() and qc.barcode.is_unique
        xy = qc[['pixel_col', 'pixel_row']].to_numpy(float)
        assert np.isfinite(xy).all()
        fig, axes = plt.subplots(2, 3, figsize=(10.8, 7.3), layout='constrained')
        rawlog = np.log10(qc.raw_UMI.to_numpy(float))
        corrlog = np.log10(qc.background_corrected_total.to_numpy(float))
        assert np.isfinite(rawlog).all() and np.isfinite(corrlog).all()
        low, high = float(min(rawlog.min(), corrlog.min())), float(max(rawlog.max(), corrlog.max()))
        for ax, ids, values, title in [(axes[0, 0], np.ones(len(qc), dtype=bool), rawlog, 'a  Original UTT: raw integer UMI'),
                                       (axes[0, 1], keep, corrlog, 'b  Primary retained: background corrected')]:
            ax.scatter(xy[:, 0], xy[:, 1], s=4, c='#dddddd', linewidths=0, rasterized=True)
            sc = ax.scatter(xy[ids, 0], xy[ids, 1], s=4, c=values[ids], vmin=low, vmax=high,
                            cmap='viridis', linewidths=0, rasterized=True)
            ax.set(title=title, xlabel='Full-resolution pixel column', ylabel='Full-resolution pixel row')
            ax.set_aspect('equal'); ax.invert_yaxis()
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_locator(MaxNLocator(5))
            ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
            ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
            cb = fig.colorbar(sc, ax=ax, shrink=.7)
            cb.set_label('log10 total (layer shown)')
        ax = axes[0, 2]
        for mask, color, name in [(~keep & ~manual, '#bcbcbc', 'QC excluded'),
                                   (manual, '#d79b5d', 'Author manual exclusion'), (keep, '#447a9b', 'Retained')]:
            ax.scatter(qc.author_Z_log10_corrected_total[mask], qc.author_Z_corrected_mito_fraction[mask],
                       s=6, c=color, label=f'{name}: {int(mask.sum())}', linewidths=0, rasterized=True)
        ax.axvline(-4, color='#333333', lw=.9); ax.axhline(4, color='#333333', lw=.9)
        ax.axvline(-3.5, color='#777777', lw=.8, ls=':'); ax.axhline(3.5, color='#777777', lw=.8, ls=':')
        ax.set(title='c  Frozen source-defined QC gates', xlabel='Z(log10 corrected total)', ylabel='Z(corrected mitochondrial fraction)')
        ax.legend(fontsize=6, loc='upper left', frameon=True, facecolor='white', framealpha=1, edgecolor='none')
        ax = axes[1, 0]
        positive = bg.background_mean.to_numpy(float)
        positive = positive[positive > 0]
        if len(positive):
            ax.hist(np.log10(positive), bins=45, color='#8da9b7')
        ax.set(title='d  ALL-minus-UTT background by original row', xlabel='log10 positive background mean', ylabel='Original gene rows')
        ax.text(.03, .97, f'{len(positive):,}/{len(bg):,} rows nonzero', transform=ax.transAxes, va='top', fontsize=7)
        ax = axes[1, 1]
        bins = np.histogram_bin_edges(np.r_[qc.raw_nGenes, qc.corrected_nGenes[keep]], bins=40)
        for values, color, label in [(qc.raw_nGenes, '#aaaaaa', 'All UTT, raw'),
                                      (qc.corrected_nGenes[keep], '#447a9b', 'Primary, corrected')]:
            ax.hist(values, bins=bins, histtype='step', color=color, lw=1.3, label=label)
        ax.set(title='e  Detected original gene rows', xlabel='Number of detected rows per spot', ylabel='Spots')
        ax.legend(fontsize=7)
        ax = axes[1, 2]
        bins = np.histogram_bin_edges(100 * np.r_[qc.raw_mito_fraction, qc.corrected_mito_fraction[keep]], bins=40)
        for values, color, label in [(100 * qc.raw_mito_fraction, '#aaaaaa', 'All UTT, raw'),
                                      (100 * qc.corrected_mito_fraction[keep], '#447a9b', 'Primary, corrected')]:
            ax.hist(values, bins=bins, histtype='step', color=color, lw=1.3, label=label)
        ax.set(title='f  Mitochondrial fraction', xlabel='Mitochondrial % (exact source name pattern)', ylabel='Spots')
        ax.legend(fontsize=7)
        fig.suptitle(f'{donor}: input UTT {len(qc):,}; primary {int(keep.sum()):,}; strict-z3.5 {int(strict.sum()):,}; manual retained 0', fontsize=10)
        for ext in ['pdf', 'png']:
            fig.savefig(out / f'{donor}_preprocessing_QC.{ext}', dpi=200, bbox_inches='tight')
        plt.close(fig)
        qc.to_csv(source_dir / f'{donor}_UTT_QC_all_spots.tsv.gz', sep='\t', index=False, compression='gzip')
        summary.append({'donor': donor, 'input_UTT_spots': len(qc), 'primary_spots': int(keep.sum()),
                        'strict_z3_5_spots': int(strict.sum()), 'manual_UTT_spots': int(manual.sum()),
                        'manual_retained': 0, 'original_gene_rows': len(bg), 'background_nonzero_gene_rows': len(positive),
                        'raw_UMI_median_all_UTT': float(qc.raw_UMI.median()),
                        'corrected_total_median_primary': float(qc.background_corrected_total[keep].median()),
                        'raw_nGenes_median_all_UTT': float(qc.raw_nGenes.median()),
                        'corrected_nGenes_median_primary': float(qc.corrected_nGenes[keep].median()),
                        'raw_mito_fraction_median_all_UTT': float(qc.raw_mito_fraction.median()),
                        'corrected_mito_fraction_median_primary': float(qc.corrected_mito_fraction[keep].median()),
                        'source_QC_path': str(paths[0].relative_to(root)), 'source_QC_sha256': sha(paths[0])})
    pd.DataFrame(summary).to_csv(out / 'DONOR_QC_DISPLAY_SUMMARY.tsv', sep='\t', index=False)
    outputs = {str(p.relative_to(root)): sha(p) for p in sorted(out.rglob('*')) if p.is_file()}
    (out / 'RENDER_STATUS.json').write_text(json.dumps({'status': 'RENDERED_PENDING_VISUAL_QA',
        'created_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'new_fit': False,
        'source_sha256': sources, 'output_sha256': outputs, 'contract': contract}, indent=2) + '\n')
    print(json.dumps({'status': 'RENDERED_PENDING_VISUAL_QA', 'donors': DONORS, 'output_directory': str(out)}), flush=True)


if __name__ == '__main__':
    main()
