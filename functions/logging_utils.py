"""
Thin helper to route verbose prints through Python's standard logging.

Call logging.basicConfig(level=logging.DEBUG) (or WARNING) in the entrypoint
to control whether these messages appear. vprint emits at DEBUG level.
"""

import logging
import inspect


def vprint(message):
    """Log message at DEBUG level using the caller's module logger."""
    try:
        frame = inspect.currentframe()
        caller = frame.f_back if frame is not None else None
        module = inspect.getmodule(caller) if caller is not None else None
        logger_name = module.__name__ if module and hasattr(module, "__name__") else __name__
        logging.getLogger(logger_name).debug(message)
    except Exception:
        logging.getLogger(__name__).debug(message)


