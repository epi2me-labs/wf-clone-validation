"""Test run_plannotate.py."""
from datetime import date, datetime
import filecmp
import os

import pytest
import workflow_glue.run_plannotate as run_plannotate


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

    # add date to template
    today = date.today()
    x = datetime.strptime(str(today), '%Y-%m-%d')
    fmt_date = x.strftime('%d-%b-%Y').upper()
    input_gbk = f"{test_data}/{gbk_file}"
    with open(input_gbk, "r") as inputfile:
        expected_gbk = f"{inputfile.read()}".format(date=fmt_date)
    with open("gbk_with_date.gbk", "w") as outfile:
        outfile.write(expected_gbk)

    # run plannotate
    run_plannotate.make_yaml('Default')
    assembly = f"{test_data}/{assembly_file}"
    run_plannotate.per_assembly(assembly, "barcode01")
    gbk = "barcode01.annotations.gbk"
    gbk_bool = filecmp.cmp(gbk, "gbk_with_date.gbk")
    os.chdir(retval)
    assert gbk_bool is True
