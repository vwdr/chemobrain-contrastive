"""Thread configuration. Import and call ``configure`` before importing numpy/torch.

Thread counts are set inside Python (not via command-line environment prefixes).
A script accepts ``--threads=N`` on its command line; the default is 3.
"""
import os
import sys


def configure(default: int = 3) -> int:
    n = default
    for arg in sys.argv[1:]:
        if arg.startswith("--threads="):
            n = int(arg.split("=", 1)[1])
    for key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[key] = str(n)
    return n


def strip_threads_arg(argv):
    return [a for a in argv if not a.startswith("--threads=")]
