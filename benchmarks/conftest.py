import pytest

# This allows you to add global configuration for your benchmarks later
def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark benchmark as slow to run")