from __future__ import annotations

import multiprocessing as _mp


_original_set_start_method = _mp.set_start_method


def _safe_set_start_method(method, force=False):
    try:
        return _original_set_start_method(method, force=force)
    except RuntimeError as exc:
        if "context has already been set" in str(exc):
            return None
        raise


_mp.set_start_method = _safe_set_start_method
