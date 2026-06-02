"""Shim — graph_path_search retired in favor of matdetangler.graph_classifier.per_k_caller.

Old import paths re-export from legacy/ for backward compatibility. New code
should use `matdetangler.run_per_k` + `matdetangler.pick_k` instead.
"""
from matdetangler.legacy.graph_path_search import *  # noqa: F401, F403
