import inspect
from unittest.mock import MagicMock, patch

import pytest

from pyquaidsce import nlsur, quaidsce
from pyquaidsce._threads import blas_thread_limit


def test_160_numerical_defaults():
    qsig = inspect.signature(quaidsce)
    nsig = inspect.signature(nlsur)
    assert qsig.parameters["sigma_tol"].default == 1e-5
    assert qsig.parameters["boot_sigma_tol"].default == 1e-5
    assert qsig.parameters["blas_threads"].default == 1
    assert nsig.parameters["sigma_tol"].default == 1e-5
    assert nsig.parameters["blas_threads"].default == 1


def test_blas_limit_is_scoped_and_none_is_unmanaged():
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


def test_blas_limit_rejects_invalid_counts():
    for bad in (0, -1, 1.5, True, "two"):
        with pytest.raises(ValueError, match="blas_threads"):
            with blas_thread_limit(bad):
                pass
