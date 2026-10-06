import inspect
import unittest
from unittest.mock import MagicMock, patch

from pyquaidsce import nlsur, quaidsce
from pyquaidsce._threads import blas_thread_limit


class ThreadAndDefaultTests(unittest.TestCase):
    def test_explicit_numerical_defaults(self):
        qsig = inspect.signature(quaidsce)
        nsig = inspect.signature(nlsur)
        self.assertEqual(qsig.parameters["outer_param_tol"].default, 1e-5)
        self.assertEqual(qsig.parameters["blas_threads"].default, 1)
        self.assertEqual(nsig.parameters["outer_param_tol"].default, 1e-5)
        self.assertEqual(nsig.parameters["blas_threads"].default, 1)

    def test_blas_limit_is_scoped_and_none_is_unmanaged(self):
        manager = MagicMock()
        manager.__enter__.return_value = None
        manager.__exit__.return_value = False

        with patch("pyquaidsce._threads.threadpool_limits", return_value=manager) as limits:
            with blas_thread_limit(1):
                pass

        limits.assert_called_once_with(limits=1, user_api="blas")
        manager.__enter__.assert_called_once()
        manager.__exit__.assert_called_once()

        with patch("pyquaidsce._threads.threadpool_limits") as limits:
            with blas_thread_limit(None):
                pass

        limits.assert_not_called()

    def test_blas_limit_rejects_invalid_counts(self):
        for bad in (0, -1, 1.5, True, "two"):
            with self.subTest(bad=bad):
                with self.assertRaisesRegex(ValueError, "blas_threads"):
                    with blas_thread_limit(bad):
                        pass


if __name__ == "__main__":
    unittest.main()
