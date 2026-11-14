PUBLIC_INTERNAL_MESSAGE = "Внутренняя ошибка сервера. Повторите запрос позже."


class ApiError(Exception):
    """Rendered as {"error": {"code": ..., "message": ...}} by main.py."""

    def __init__(self, status: int, code: str, message: str):
        super().__init__(message)
        self.status, self.code, self.message = status, code, message
