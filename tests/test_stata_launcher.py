"""Regression tests for the Stata-embedded Python subprocess launcher.

These use a fake Stata SFI module; no Stata license or long estimation is needed.
"""

from __future__ import annotations

import os
import sys
import tempfile
import types
import unittest
from unittest import mock

import pyquaidsce
from pyquaidsce import stata_bridge


class StataLauncherTests(unittest.TestCase):
    def test_get_executable_uses_stata_setting(self):
        sfi = types.SimpleNamespace(
            Macro=types.SimpleNamespace(
                getGlobal=lambda name: '"/opt/python with space/bin/python"'
                if name == "c(python_exec)" else ""
            )
        )
        with mock.patch.dict(sys.modules, {"sfi": sfi}), mock.patch(
            "os.path.isfile", return_value=True
        ):
            self.assertEqual(
                stata_bridge._get_stata_python_executable(),
                "/opt/python with space/bin/python",
            )

    def test_missing_executable_is_rejected(self):
        sfi = types.SimpleNamespace(
            Macro=types.SimpleNamespace(getGlobal=lambda name: "")
        )
        with mock.patch.dict(sys.modules, {"sfi": sfi}):
            with self.assertRaisesRegex(RuntimeError, "no valid Python"):
                stata_bridge._get_stata_python_executable()

    def test_child_import_identity_is_checked(self):
        path = os.path.realpath(pyquaidsce.__file__)
        completed = types.SimpleNamespace(
            returncode=0,
            stdout=f"{path}\n{pyquaidsce.__version__}\n",
            stderr="",
        )
        with mock.patch("subprocess.run", return_value=completed) as run:
            stata_bridge._verify_child_python("/tmp/py")
            self.assertEqual(run.call_args.args[0][0], "/tmp/py")
            self.assertIn("-c", run.call_args.args[0])

        wrong = types.SimpleNamespace(
            returncode=0, stdout=f"{path}\n0.0.0\n", stderr=""
        )
        with mock.patch("subprocess.run", return_value=wrong):
            with self.assertRaisesRegex(RuntimeError, "different pyquaidsce"):
                stata_bridge._verify_child_python("/tmp/py")

        missing = types.SimpleNamespace(
            returncode=1, stdout="", stderr="No module named pyquaidsce"
        )
        with mock.patch("subprocess.run", return_value=missing):
            with self.assertRaisesRegex(RuntimeError, "cannot import pyquaidsce"):
                stata_bridge._verify_child_python("/tmp/py")

    def test_launch_uses_selected_python_not_sys_executable(self):
        with tempfile.TemporaryDirectory() as temp, mock.patch(
            "tempfile.gettempdir", return_value=temp
        ), mock.patch.object(
            stata_bridge, "_get_stata_python_executable", return_value="/chosen/python"
        ), mock.patch.object(
            stata_bridge, "_verify_child_python"
        ), mock.patch(
            "subprocess.Popen"
        ) as popen:
            try:
                stata_bridge._launch_job({"data": [1]}, {"reps": 0})
                command = popen.call_args.args[0]
                self.assertEqual(command[0], "/chosen/python")
                self.assertEqual(command[1:3], ["-m", "pyquaidsce.bootstrap_runner"])
                self.assertNotEqual(command[0], sys.executable)
                self.assertIn("stderr_path", stata_bridge._BOOT_STATE)
            finally:
                stata_bridge._cleanup_boot_files()

    def test_exited_child_surfaces_its_stderr(self):
        messages = {}
        sfi = types.SimpleNamespace(
            Scalar=types.SimpleNamespace(setValue=lambda key, value: messages.__setitem__(key, value)),
            Macro=types.SimpleNamespace(setLocal=lambda key, value: messages.__setitem__(key, value)),
        )
        with tempfile.TemporaryDirectory() as temp:
            stderr_file = os.path.join(temp, "stderr.txt")
            with open(stderr_file, "w", encoding="utf-8") as output:
                output.write("No module named pyquaidsce.bootstrap_runner\n")
            stata_bridge._BOOT_STATE.update(
                result_path=os.path.join(temp, "absent_result.pkl"),
                status_path="",
                stderr_path=stderr_file,
                job_path="",
                proc=types.SimpleNamespace(poll=lambda: 1, returncode=1),
            )
            try:
                with mock.patch.dict(sys.modules, {"sfi": sfi}):
                    stata_bridge.poll_bootstrap()
                self.assertEqual(messages["_pyq_boot_err"], 1)
                self.assertIn(
                    "No module named pyquaidsce.bootstrap_runner",
                    messages["_pyq_boot_errmsg"],
                )
            finally:
                stata_bridge._BOOT_STATE.clear()


if __name__ == "__main__":
    unittest.main()
