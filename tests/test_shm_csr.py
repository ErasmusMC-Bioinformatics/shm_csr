# Copyright (c) 2021 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

GIT_ROOT = str(Path(__file__).parent.parent.absolute())
TEST_DIR = Path(__file__).parent
TEST_DATA_DIR = TEST_DIR / "data"
VALIDATION_DATA_DIR = TEST_DIR / "validation_data"
CONTROL_NWK377_PB_IGHC_MID1_40nt_2 = TEST_DATA_DIR / "CONTROL_NWK377_PB_IGHC_MID1_40nt_2.txz"


@pytest.fixture(scope="module")
def shm_csr_result():
    temp_dir = tempfile.mktemp()
    shutil.copytree(GIT_ROOT, temp_dir)
    input = str(CONTROL_NWK377_PB_IGHC_MID1_40nt_2)
    out_files_path = os.path.join(temp_dir, "results")
    out_file = os.path.join(out_files_path, "result.html")
    infile_name = "input_data"
    functionality = "productive"
    unique = "Sequence.ID"
    naive_output = "no"
    naive_output_ca = "None"
    naive_output_cg = "None"
    naive_output_cm = "None"
    naive_output_ce = "None"
    naive_output_all = "None"
    filter_unique = "remove"
    filter_unique_count = '2'
    class_filter = '70_70'
    empty_region_filter = 'FR1'
    fast = 'no'
    cmd = [
        "bash",
        "wrapper.sh",
        input,
        "custom",
        out_file,
        out_files_path,
        infile_name,
        "-",
        functionality,
        unique,
        naive_output,
        naive_output_ca,
        naive_output_cg,
        naive_output_cm,
        naive_output_ce,
        naive_output_all,
        filter_unique,
        filter_unique_count,
        class_filter,
        empty_region_filter,
        fast
    ]
    subprocess.run(cmd, cwd=temp_dir, stdout=sys.stdout, stderr=sys.stderr,
                   check=True)
    yield Path(out_files_path)
    shutil.rmtree(temp_dir)


def test_check_output(shm_csr_result):
    assert shm_csr_result.exists()


@pytest.mark.parametrize("filename", [
    "aa_histogram_sum.txt",
])
def test_aa_histogram_sum(shm_csr_result, filename):
    with open(Path(shm_csr_result, filename)) as result_h:
        with open(Path(VALIDATION_DATA_DIR, filename)) as validate_h:
            for line in result_h:
                assert line == validate_h.readline()
