"""Wrapper CLI lint: every flag the bash wrapper passes to a matdetangler.<mod>
python module must be a flag that module's argparse accepts.

Catches the class of bug where a wrapper edit forwards a flag to the wrong
module (e.g. --include-flanks on graph_path_search.py, --cache-dir on a module
that has dropped that flag, etc.). Unit tests on individual python modules
test each module in isolation — they NEVER see the wrapper's flag-routing
choices. This test bridges the gap.

Detection method:
  1. Read the bash wrapper file.
  2. Find every `python -m matdetangler.<mod> ...` invocation. The invocation
     spans multiple bash lines because it uses backslash continuations.
  3. Extract the set of long flags (`--foo`) appearing in that span.
  4. Import the python module and inspect its argparse parser; collect the
     set of long flags it accepts.
  5. Assert (flags_in_wrapper - flags_accepted_by_module) == set().

Limitations:
  - Bash conditionals like `$( [ "$DEBUG" -eq 1 ] && echo --foo )` are still
    flags that COULD be forwarded; we lint them as if always present
    (correct behavior — they need to be accepted).
  - Short flags (`-t`) and positional args are not linted (the wrapper uses
    only long flags for python module invocations).
  - The wrapper's own flags (those the wrapper parses for itself before
    dispatching to a module) are not linted by this test.
"""
import os, re, sys, unittest
import argparse

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

WRAPPER = os.path.join(ROOT, "MATdetangler-cli")


def _read_wrapper() -> str:
    with open(WRAPPER) as f:
        return f.read()


# Join multi-line python-invocation blocks. Bash continues lines with a trailing
# backslash; we collapse those into single logical lines.
def _join_continuations(src: str) -> str:
    return re.sub(r"\\\s*\n\s*", " ", src)


# Match every `python -m matdetangler.<mod>  ... ` block up to (but not
# including) the next pipe `|`, redirect `>`, semicolon `;`, or end-of-line.
_INVOKE_RE = re.compile(
    r'(?:"\$PY"|python3?(?:\.\d+)?)\s+-m\s+matdetangler\.(\w+)\s+([^\n|;>]*)')
_LONG_FLAG_RE = re.compile(r'--([\w-]+)')


def _flags_passed_to(src_joined: str, mod_name: str) -> set[str]:
    """Return all long flag names (without leading --) that the wrapper passes
    to python module `mod_name` across all its invocations."""
    flags: set[str] = set()
    for m in _INVOKE_RE.finditer(src_joined):
        if m.group(1) != mod_name: continue
        flags.update(_LONG_FLAG_RE.findall(m.group(2)))
    return flags


def _argparse_long_options(module_name: str) -> set[str]:
    """Import matdetangler.<module_name>, build its CLI parser, return the set
    of long option strings the parser accepts (without the leading '--').

    The convention in this codebase is: each CLI module defines a `_cli(argv)`
    function that builds a parser and parses argv. We can't safely call _cli
    (it has side effects), so we monkey-patch parse_args to capture the
    parser without parsing.
    """
    import importlib
    mod = importlib.import_module(f"matdetangler.{module_name}")
    if not hasattr(mod, "_cli"):
        raise RuntimeError(f"matdetangler.{module_name} has no _cli() entry")

    captured: dict = {}
    orig_parse_args = argparse.ArgumentParser.parse_args

    def fake_parse_args(self, args=None, namespace=None):
        captured["parser"] = self
        # Return a Namespace with every dest set to None — enough for callers
        # that immediately dispatch (we won't let them, see below).
        ns = argparse.Namespace()
        for a in self._actions:
            if a.dest != "help":
                setattr(ns, a.dest, None)
        raise _CapturedParser()    # abort _cli before it does anything

    class _CapturedParser(Exception): pass

    argparse.ArgumentParser.parse_args = fake_parse_args
    try:
        try:
            mod._cli([])
        except _CapturedParser:
            pass
        except SystemExit:
            pass    # argparse may sys.exit on bad args; we caught the parser first
    finally:
        argparse.ArgumentParser.parse_args = orig_parse_args

    parser = captured.get("parser")
    if parser is None:
        raise RuntimeError(f"failed to capture argparse parser from matdetangler.{module_name}._cli")

    def _collect(p, into: set[str]) -> None:
        for action in p._actions:
            for opt in action.option_strings:
                if opt.startswith("--"):
                    into.add(opt[2:])
            # Descend into subparsers (e.g. matdetangler.cluster has `align`/`cut`).
            if isinstance(action, argparse._SubParsersAction):
                for subp in action.choices.values():
                    _collect(subp, into)

    long_opts: set[str] = set()
    _collect(parser, long_opts)
    return long_opts


class TestWrapperFlags(unittest.TestCase):
    """One test per python module the wrapper invokes."""

    @classmethod
    def setUpClass(cls):
        cls.src = _join_continuations(_read_wrapper())

    def _check(self, mod_name: str) -> None:
        passed = _flags_passed_to(self.src, mod_name)
        if not passed:
            self.skipTest(f"wrapper never invokes matdetangler.{mod_name}")
        accepted = _argparse_long_options(mod_name)
        unknown = passed - accepted
        self.assertFalse(
            unknown,
            f"Wrapper passes flag(s) to matdetangler.{mod_name} that its argparse rejects: "
            f"{sorted(unknown)}.  Accepted flags: {sorted(accepted)}")

    def test_input_process_flags(self):    self._check("input_process")
    def test_anchor_search_flags(self):    self._check("anchor_search")
    def test_graph_path_search_flags(self):self._check("graph_path_search")
    def test_pick_alleles_flags(self):     self._check("pick_alleles")
    def test_graph_paths_flags(self):      self._check("graph_paths")
    def test_pairwise_identity_flags(self):self._check("pairwise_identity")
    def test_consensus_qc_flags(self):     self._check("consensus_qc")
    def test_summary_table_flags(self):    self._check("summary_table")
    def test_cluster_flags(self):          self._check("cluster")


class TestWrapperLintSanity(unittest.TestCase):
    """Sanity checks on the lint mechanism itself — make sure the regex and
    argparse-capture actually work on known examples."""

    def test_finds_at_least_one_invocation(self):
        src = _join_continuations(_read_wrapper())
        n_invocations = len(_INVOKE_RE.findall(src))
        self.assertGreater(n_invocations, 5,
                           "wrapper invokes matdetangler modules but the regex found <=5 — broken matcher")

    def test_captures_a_known_flag(self):
        """anchor_search definitely has --sample; assert the capture sees it."""
        opts = _argparse_long_options("anchor_search")
        self.assertIn("sample", opts)
        self.assertIn("source", opts)
        self.assertIn("spades-dir", opts)

    def test_long_flag_re_extracts_basic(self):
        s = "  --sample foo --threads 8 --include-flanks --blast-out-dir bar"
        self.assertEqual(set(_LONG_FLAG_RE.findall(s)),
                         {"sample", "threads", "include-flanks", "blast-out-dir"})


if __name__ == "__main__":
    unittest.main(verbosity=2)
