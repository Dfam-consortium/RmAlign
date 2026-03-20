"""
telib/rmlog.py — Centralised logging setup for RmAlign pipeline tools.

Call ``setup(level_str, log_file)`` once at the start of each script's
``main()``.  All modules then use the standard idiom::

    import logging
    logger = logging.getLogger(__name__)
    logger.debug("...")

Log routing
-----------
* **With a log file** (Nextflow mode):
  - File handler captures at *level_str* and above.
  - stderr handler captures WARNING and above only, so INFO/DEBUG output
    does not pollute Nextflow's ``.command.err`` (which RepeatMasker scans
    after a failed run to surface errors).

* **Without a log file** (standalone / interactive mode):
  - stderr handler captures at *level_str*, so debug output is visible
    directly.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional

_FMT = "%(asctime)s [%(levelname)-8s] %(name)s: %(message)s"
_DATEFMT = "%Y-%m-%dT%H:%M:%S"

_LEVELS = {"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"}


def setup(level_str: str = "INFO", log_file: Optional[str] = None) -> None:
    """Configure the root logger for an RmAlign pipeline script.

    Parameters
    ----------
    level_str:
        One of DEBUG / INFO / WARNING / ERROR / CRITICAL (case-insensitive).
        Unknown strings fall back to INFO.
    log_file:
        Optional path for the per-task log file.  The parent directory is
        created if it does not exist.
    """
    level = getattr(logging, level_str.upper(), logging.INFO)
    fmt = logging.Formatter(_FMT, datefmt=_DATEFMT)

    handlers: list[logging.Handler] = []

    if log_file:
        # --- Nextflow / file mode ---
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_file, mode="a", encoding="utf-8", delay=True)
        fh.setLevel(level)
        fh.setFormatter(fmt)
        handlers.append(fh)

        # Keep stderr quiet (WARNING+) so it stays clean for error detection
        sh = logging.StreamHandler(sys.stderr)
        sh.setLevel(logging.WARNING)
        sh.setFormatter(fmt)
        handlers.append(sh)
    else:
        # --- Standalone / interactive mode ---
        sh = logging.StreamHandler(sys.stderr)
        sh.setLevel(level)
        sh.setFormatter(fmt)
        handlers.append(sh)

    root = logging.getLogger()
    root.setLevel(level)
    for h in handlers:
        root.addHandler(h)


def level_from_str(s: str) -> int:
    """Return the ``logging`` int for a level name string, defaulting to INFO."""
    return getattr(logging, s.upper(), logging.INFO)
