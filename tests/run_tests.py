# Copyright 2026 Project Contributors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Standalone CLI test runner executing pytest and returning clean exit codes."""

from __future__ import annotations

import os
from pathlib import Path
import sys

os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
sys.dont_write_bytecode = True

import pytest


def main() -> int:
    """Runs the test suite in tests/test_jaspar_visualize.py."""
    test_file = Path(__file__).resolve().parent / "test_jaspar_visualize.py"
    return int(pytest.main(["-v", "-p", "no:cacheprovider", str(test_file)]))


if __name__ == "__main__":
    sys.exit(main())
