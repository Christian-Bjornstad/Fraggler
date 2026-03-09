"""
Runtime hook: monkey-patches Python 3.10's dis._get_const_info to handle
the IndexError that occurs when scanning certain scipy/sklearn .pyc files.

This is a workaround for a known Python 3.10.0 bug in the dis module.
"""
import dis

_original_get_const_info = dis._get_const_info


def _patched_get_const_info(const_index, const_list):
    try:
        return _original_get_const_info(const_index, const_list)
    except IndexError:
        return const_index, repr(const_index)


dis._get_const_info = _patched_get_const_info
