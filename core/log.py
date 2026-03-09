"""
Fraggler Diagnostics — Centralized Logging

Provides a param-based LogBuffer that the GUI can watch, and a global
log() function for all modules.
"""
from __future__ import annotations

import param
import logging
from datetime import datetime

logger = logging.getLogger()

class LogBuffer(param.Parameterized):
    """Observable text buffer for the GUI log viewer."""
    text = param.String(default="")

    def write(self, msg: str) -> None:
        self.text += str(msg) + "\n"
        self.param.trigger("text")

    def clear(self) -> None:
        self.text = ""
        self.param.trigger("text")


# Singleton
log_buffer = LogBuffer()


def log(msg: str) -> None:
    """Append a timestamped message to the log buffer and stdout."""
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    log_buffer.write(line)
    logger.info(msg)
