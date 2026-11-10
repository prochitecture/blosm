"""Select runtime code without changing the location of test sources/fixtures."""
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ADDON_ROOT = Path(os.environ.get('BLOSM_TEST_ROOT', REPO)).resolve()
PACKAGE = os.environ.get('BLOSM_TEST_PACKAGE', REPO.name)
PROFILE = os.environ.get('BLOSM_TEST_PROFILE', 'general')


def assert_target(addon):
    if Path(addon.__file__).resolve().parent != ADDON_ROOT:
        raise AssertionError('Test loaded an addon outside the selected candidate')
    if PROFILE != 'general' and ADDON_ROOT.is_relative_to(REPO):
        raise AssertionError('A build profile must not load the development checkout')
