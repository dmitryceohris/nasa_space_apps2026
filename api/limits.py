"""Single-process, bounded per-peer token buckets for expensive API endpoints."""
from collections import OrderedDict
import math
import threading
import time


class RequestLimiter:
    def __init__(self, per_minute=120, max_clients=4096, clock=None):
        if per_minute < 1 or max_clients < 1:
            raise ValueError('Request limits must be positive')
        self.capacity = per_minute
        self.refill = per_minute / 60
        self.max_clients = max_clients
        self.clock = clock or time.monotonic
        self.clients = OrderedDict()
        self.lock = threading.Lock()

    def retry_after(self, peer):
        with self.lock:
            now = self.clock()
            tokens, previous = self.clients.pop(peer, (self.capacity, now))
            tokens = min(self.capacity, tokens + max(0, now - previous) * self.refill)
            delay = max(1, math.ceil((1 - tokens) / self.refill)) if tokens < 1 else None
            self.clients[peer] = (tokens if delay else tokens - 1, now)
            while len(self.clients) > self.max_clients:
                self.clients.popitem(last=False)
            return delay
