"""Test find_inserts.py."""
import os

import pandas as pd
import pysam
import pytest
import workflow_glue.find_inserts as find_inserts

# Use FASTA files as expected inserts.
# The same sequences used for gitlab CI test.
EXPECTED_INSERTS = [
    "barcode01.insert.fasta",
    "barcode02.insert.fasta",
    "barcode03.insert.fasta",
    "barcode04.insert.fasta"]


def get_fasta_seq_dic(file_path, file_list):
    """Get list of fasta sequences as dictionary."""
    fasta_dic = {}
    for file_n in file_list:
        file_name = f"{file_path}/{file_n}"
        with pysam.FastaFile(file_name) as f:
            barcode = f.references[0]
            seq = f.fetch(barcode)
            fasta_dic[barcode] = seq
    return fasta_dic


@pytest.fixture
def test_data(request):
    """Define data location fixture."""
    return os.path.join(
        request.config.getoption("--test_data"),
        "workflow_glue",
        "find_inserts")


@pytest.mark.parametrize(
        "test_data_fixture, bed_files, assemblies_dir, \
            expected_inserts, expected_df_file",
        [(
            "test_data", "insert_beds",
            "assemblies", EXPECTED_INSERTS, "expected_df.csv")]
        )
def test_read_seqkit(
        test_data_fixture, bed_files, assemblies_dir,
        expected_inserts, expected_df_file, request):
    """Test insert sequence generation from seqkit amplicon output."""
    # Prepare expected insert dictionary and data frame
    test_data_dir = request.getfixturevalue(test_data_fixture)
    bed_files = f"{test_data_dir}/{bed_files}"
    expected_fp = f"{test_data_dir}/expected_insert"
    expected_seqs_dic = get_fasta_seq_dic(expected_fp, expected_inserts)
    exp_df_file = f"{test_data_dir}/{expected_df_file}"
    expected_df = pd.read_csv(exp_df_file).sort_values('Sample')

    # Get actual insert dictionary and data frame
    assemblies_dir = f"{test_data_dir}/{assemblies_dir}"
    sk_amplicon_beds = [
        os.path.join(bed_files, file) for file in os.listdir(bed_files)]
    actual_df, actual_seqs_dic = find_inserts.read_seqkit(
        sk_amplicon_beds, assemblies_dir=assemblies_dir)

    # Compare
    pd.testing.assert_frame_equal(
        expected_df,
        actual_df.sort_values('Sample').reset_index(drop=True),
        check_index_type=False)
    for k, v in expected_seqs_dic.items():
        assert actual_seqs_dic[k] == v


@pytest.mark.parametrize(
        "seq, expected",
        [
            ("GGGATATAGCCCCGCATAT", "ATATGCGGGGCTATATCCC"),
            ("TATCCCGCCCCCXCAGCTTGCCAGNTCTTT",
             "AAAGANCTGGCAAGCTGXGGGGGCGGGATA")
        ]
        )
def test_reverse_complement(seq, expected):
    """Test reverse complement."""
    actual = find_inserts.reverse_complement(seq)
    assert actual == expected
