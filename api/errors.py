class ApiError(Exception):
    """Rendered as {"error": {"code": ..., "message": ...}} by main.py."""

    def __init__(self, status: int, code: str, message: str, headers: dict[str, str] | None = None):
        super().__init__(message)
        self.status, self.code, self.message, self.headers = status, code, message, headers
