"""report tests."""
import os

import pandas as pd
import pytest
from workflow_glue.report_utils.report_utils import get_cutsite_table


@pytest.fixture
def test_data(request):
    """Define data location fixture."""
    return os.path.join(
        request.config.getoption("--test_data"),
        "workflow_glue",
        "report")


@pytest.mark.parametrize(
        "cutsite_csv,samples,expected_values",
        [
            ("cut_sites.csv",
             ["sample01", "sample02", "sample03", "sample04"],
             [13.61, 28.86, 75.0, "N/A"])
        ]
    )
def test_get_cutsite_table(test_data, cutsite_csv, samples, expected_values):
    """Test get cutsite table."""
    expected_df = pd.DataFrame(
        {"Sample": samples, "Linearisation efficiency (%)": expected_values}
        )
    cutsite_file = f"{test_data}/{cutsite_csv}"
    actual_df = get_cutsite_table(cutsite_file, samples)
    pd.testing.assert_frame_equal(
        actual_df.reset_index(drop=True), expected_df.reset_index(drop=True))
