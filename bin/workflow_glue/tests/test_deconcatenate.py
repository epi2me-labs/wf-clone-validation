"""Test deconcatenate.py."""

import deconcatenate
import mappy as mp
import pytest


@pytest.mark.parametrize("test_data, approx_size, expected", [
    # is approx size so don't remove repeat
    ("test_data/workflow_glue/deconcatenate/barcode01.fasta",
        4690, "test_data/workflow_glue/deconcatenate/barcode01_expected.fasta"),
    # single repeat
    ("test_data/workflow_glue/deconcatenate/barcode02.fasta",
        2345, "test_data/workflow_glue/deconcatenate/barcode02_expected.fasta"),
    ("test_data/workflow_glue/deconcatenate/barcode03.fasta",
        2345, "test_data/workflow_glue/deconcatenate/barcode03_expected.fasta"),
    # multiple repeats, remove one
    ("test_data/workflow_glue/deconcatenate/barcode04.fasta",
        4690, "test_data/workflow_glue/deconcatenate/barcode04_expected.fasta"),
    # no self alignment
    ("test_data/workflow_glue/deconcatenate/barcode05.fasta",
        1000, "test_data/workflow_glue/deconcatenate/barcode05_expected.fasta")
])
def test_deconcatenate(test_data, approx_size, expected):
    """Test deconcatenate."""
    test_seq_list = mp.fastx_read(test_data)
    expected_seq_list = mp.fastx_read(expected)
    for test_seq, exp_seq in zip(test_seq_list, expected_seq_list):
        corrected = deconcatenate.deconcatenate(test_seq[1], approx_size)
        assert corrected == exp_seq[1]
