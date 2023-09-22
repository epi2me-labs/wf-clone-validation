"""Test run_plannotate.py."""
import filecmp
import os

import pytest
from workflow_glue import run_plannotate


EXPECTED = [
    ("barcode01.fasta", "barcode01.annotations.gbk")
]


@pytest.fixture
def test_data(request):
    """Define data location fixture."""
    return os.path.join(
        request.config.getoption("--test_data"),
        "workflow_glue",
        "run_plannotate")


@pytest.mark.parametrize("assembly_file,gbk_file", EXPECTED)
def test_per_assembly(test_data, assembly_file, gbk_file, tmpdir):
    """Test per assembly function in run_plannotate outputs correct gbk."""
    retval = os.getcwd()
    os.chdir(tmpdir)
    run_plannotate.make_yaml('Default')
    assembly = f"{test_data}/{assembly_file}"
    run_plannotate.per_assembly(assembly, "barcode01")
    expected_gbk = f"{test_data}/{gbk_file}"
    gbk = f"{tmpdir}/barcode01.annotations.gbk"
    os.chdir(retval)
    assert filecmp.cmp(gbk, expected_gbk) is True
