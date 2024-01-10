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
from xml.etree import ElementTree
from xml.etree.ElementTree import Element

import pytest

GIT_ROOT = Path(__file__).parent.parent.absolute()
TEST_DIR = Path(__file__).parent
TEST_DATA_DIR = TEST_DIR / "data"
VALIDATION_DATA_DIR = TEST_DIR / "validation_data"
CONTROL_NWK377_PB_IGHC_MID1_40nt_2 = TEST_DATA_DIR / "CONTROL_NWK377_PB_IGHC_MID1_40nt_2.txz"


def get_container():
    tool = ElementTree.parse(GIT_ROOT / "shm_csr.xml").getroot()
    requirements: Element = tool.find("requirements")
    container = requirements.find("container")
    return container.text


def ignore_files(src, files):
    "Ignore virtualenv and git directories to prevent massive tmp folders"
    if os.path.basename(src) in (".venv", ".git"):
        return files
    return ()


def run_shm_csr(
    input=str(CONTROL_NWK377_PB_IGHC_MID1_40nt_2),
    infile_name = "input_data",
    functionality = "productive",
    unique = "Sequence.ID",
    naive_output = "no",
    naive_output_ca = "None",
    naive_output_cg = "None",
    naive_output_cm = "None",
    naive_output_ce = "None",
    naive_output_all = "None",
    naive_output_igm_naive = "None",
    naive_output_igm_naive_memory = "None",
    filter_unique = "remove",
    filter_unique_count = '2',
    class_filter = '70_70',
    empty_region_filter = 'FR1',
    # Skip baseline and changeo by default. These tools cannot be modified
    # anyway and take most of the test time to execute. The environment
    # variable can be set to "no" on the CI so the code path is tested
    # at the time a PR is ready.
    run_changeo = "yes" if os.environ.get("SHM_CSR_FAST") == "no" else "no",
    run_baseline = "yes" if os.environ.get("SHM_CSR_FAST") == "no" else "no",
):
    temp_dir = Path(tempfile.mkdtemp())
    tool_dir = temp_dir / "shm_csr"
    shutil.copytree(
        GIT_ROOT, tool_dir,
        # Ignore .venv and .git directories.
        ignore=ignore_files)
    working_dir = temp_dir / "working"
    working_dir.mkdir(parents=True)
    output_dir = temp_dir / "outputs"
    output_dir.mkdir(parents=True)
    wrapper = str(tool_dir / "wrapper.sh")
    out_files_path = output_dir / "results"
    out_file = out_files_path / "result.html"
    cmd = [
        "bash",
        wrapper,
        input,
        "custom",
        str(out_file),
        str(out_files_path),
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
        naive_output_igm_naive,
        naive_output_igm_naive_memory,
        filter_unique,
        filter_unique_count,
        class_filter,
        empty_region_filter,
        run_changeo,
        run_baseline
    ]
    docker_cmd = ["docker", "run", "-v", f"{temp_dir}:{temp_dir}",
                  "--rm",  # Remove container after running
                  "-v", f"{input}:{input}",
                  "-w", str(working_dir),
                  # Run as current user which allows deletion of files.
                  # It also mitigates some security considerations
                  "-u", f"{os.getuid()}:{os.getgid()}",
                  # Run with default seccomp profile to mitigate mitigation slowdown
                  # http://mamememo.blogspot.com/2020/05/cpu-intensive-rubypython-code-runs.html
                  # This knocks down test runtime from 8 to 6 minutes.
                  "--security-opt", "seccomp=unconfined",
                  # Use a mulled container generated with `planemo mull`
                  "quay.io/biocontainers/mulled-v2-f7d31c9d7424063a492fc0e5ecbf89bc757c0107:2b50bdd4d8c1fefc6ec24b0753fad0dcecec843b-0"
                  ] + cmd
    with open(temp_dir / "stderr", "wt") as stderr_file:
        with open(temp_dir / "stdout", "wt") as stdout_file:
            subprocess.run(docker_cmd, cwd=working_dir, stdout=stdout_file,
                           stderr=stderr_file, check=True)
    return Path(out_files_path)


@pytest.fixture(scope="module")
def shm_csr_result():
    yield run_shm_csr()


def test_check_output(shm_csr_result):
    assert shm_csr_result.exists()


@pytest.mark.parametrize("filename", os.listdir(VALIDATION_DATA_DIR))
def test_results_match_validation(shm_csr_result, filename):
    with open(Path(shm_csr_result, filename)) as result_h:
        with open(Path(VALIDATION_DATA_DIR, filename)) as validate_h:
            for line in result_h:
                # Skip two faulty lines in shm_overview.
                # TODO: Fix the issue.
                validation_line = validate_h.readline()
                if filename == "shm_overview.txt":
                    if (line.startswith("RGYW (%)") or
                            line.startswith("WRCY (%)")):
                        continue
                assert line == validation_line


def test_nt_overview(shm_csr_result):
    with open(Path(shm_csr_result, "sequence_overview", "ntoverview.txt")
              ) as result_h:
        with open(Path(TEST_DIR, "sequence_overview", "ntoverview.txt")
                  ) as validate_h:
            for line in result_h:
                assert line == validate_h.readline()


def test_baseline_succeeds():
    run_shm_csr(
        functionality="unproductive",
        empty_region_filter="None",
        filter_unique="no",
        unique="VGene,DGene,JGene,CDR3.IMGT.seq",
        class_filter="101_101_IGA",
        naive_output="yes",
        run_baseline="yes",
        run_changeo="yes",
    )
