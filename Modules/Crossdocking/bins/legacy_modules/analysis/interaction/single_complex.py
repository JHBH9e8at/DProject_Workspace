"""Standard compatibility boundary for shared single-complex analysis."""

from ..interaction_common.single_complex import main as _shared_main
from ..interaction_common.single_complex import run_single_complex as _shared_run


def run_single_complex(*args, **kwargs):
    kwargs["_include_elapsed"] = False
    return _shared_run(*args, **kwargs)


def main(argv=None):
    return _shared_main(argv, include_elapsed=False)


if __name__ == "__main__":
    main()
