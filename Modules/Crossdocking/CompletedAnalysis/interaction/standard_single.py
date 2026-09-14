"""Sequential single-complex compatibility adapter."""

from .common.single_complex import main as _main
from .common.single_complex import run_single_complex as _run


def run_single_complex(*args, **kwargs):
    kwargs["_include_elapsed"] = False
    return _run(*args, **kwargs)


def main(argv=None):
    return _main(argv, include_elapsed=False)
