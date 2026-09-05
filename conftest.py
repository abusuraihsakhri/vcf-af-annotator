"""Pytest configuration and shared fixtures."""
import os
import pytest


@pytest.fixture(autouse=True)
def _set_audit_key(monkeypatch):
    """Set a deterministic AUDIT_SECRET_KEY for tests that exercise the audit trail."""
    monkeypatch.setenv("AUDIT_SECRET_KEY", "test-audit-key-do-not-use-in-production")
    yield
