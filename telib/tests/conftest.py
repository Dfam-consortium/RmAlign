from pathlib import Path
import pytest

# Fixture to initialize the location of test data
# files for use in tests.  E.g.
#
#   def test_something(test_data_dir):
#       foo = test_data_dir / "foo.txt"
#
@pytest.fixture(scope="session")
def test_data_dir() -> Path:
    # tests/data
    return Path(__file__).resolve().parent / "data"
