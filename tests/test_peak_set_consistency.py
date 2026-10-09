"""
The "filtered" peak path recorded in the results DB and the DB's stat columns
must describe the peak set MEME/FIMO actually consume.
"""

import report as rp
import update_db as m


def test_meme_peaks_path_without_extra_filters_is_fold_file():
    assert m.meme_peaks_path("/out/run/", "TF1", 2, "") == "/out/run/MACS/TF1_peaks_fold2.narrowPeak"


def test_meme_peaks_path_includes_enabled_filter_suffix():
    assert (m.meme_peaks_path("/out/run", "TF1", 3, "_bl_rmsk")
            == "/out/run/MACS/TF1_peaks_fold3_bl_rmsk.narrowPeak")


def test_db_stores_filter_cascade_and_frip_columns():
    for col in ("num_peaks_bl", "num_peaks_rmsk", "frip", "frip_filt"):
        assert col in m.COLS
        assert col in rp.make_cols()   # every stored column is produced by the report


def test_db_has_no_column_the_report_never_produces():
    produced = set(rp.make_cols())
    run_level = {
        "run_date", "output_dir", "genome_ref", "genome_size", "control", "threads",
        "mapq", "max_frags", "macs3_format", "macs3_foldch_levels",
        "macs3_meme_foldch_level", "meme_nmotifs", "meme_minw", "meme_maxw",
        "meme_maxpeaks", "fimo_thresh", "sample", "r1", "r2", "is_treatment",
    }
    assert set(m.COLS) - run_level <= produced


def test_report_header_says_which_columns_use_the_meme_set():
    header = rp._report_header_html(filter_foldch=5)
    assert "num_peaks_filt" in header and "5" in header
    assert "final set fed to MEME/FIMO" in header
