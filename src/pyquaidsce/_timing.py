"""Private cooperative-deadline helpers used by bootstrap replications."""

from __future__ import annotations

import time
from typing import Optional


def check_deadline(deadline: Optional[float]) -> None:
    """Raise ``TimeoutError`` after a monotonic-clock deadline.

    The estimator checks at numerical iteration and chunk boundaries.  This is
    deliberately cooperative: it works in sequential and multiprocessing
    bootstrap modes without signals or unsafe process termination.
    """
    if deadline is not None and time.perf_counter() > deadline:
        raise TimeoutError("bootstrap replication exceeded its time limit")
