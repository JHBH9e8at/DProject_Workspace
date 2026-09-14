"""Compatibility wrapper for :mod:`run_crossdocking`."""

try:
    from .run_crossdocking import *  # noqa: F401,F403
    from .run_crossdocking import main
except ImportError:
    from run_crossdocking import *  # type: ignore  # noqa: F401,F403
    from run_crossdocking import main  # type: ignore


if __name__ == "__main__":
    raise SystemExit(main())
