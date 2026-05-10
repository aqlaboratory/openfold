# Copyright 2021 AlQuraishi Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import py_compile
import unittest
import warnings


class TestScriptUtils(unittest.TestCase):
    def test_script_utils_compiles_without_syntax_warnings(self):
        repo_root = os.path.dirname(os.path.dirname(__file__))
        script_utils_path = os.path.join(
            repo_root,
            "openfold",
            "utils",
            "script_utils.py",
        )

        with warnings.catch_warnings():
            warnings.simplefilter("error", SyntaxWarning)
            py_compile.compile(script_utils_path, doraise=True)
