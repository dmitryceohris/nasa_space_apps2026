"""Keep short-interval offline polling tests independent of production request budgets."""
import sys
import pytest


@pytest.fixture(autouse=True)
def api_request_budget(monkeypatch):
    main = sys.modules.get('api.main')
    if main is not None:
        from api.limits import RequestLimiter
        monkeypatch.setattr(main, 'request_limiter', RequestLimiter(10000))
