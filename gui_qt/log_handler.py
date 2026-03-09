import logging
from PyQt6.QtCore import QObject, pyqtSignal

class QLogSignalEmitter(QObject):
    log_signal = pyqtSignal(str)

class QLogHandler(logging.Handler):
    """
    Custom logging handler that emits log messages via a PyQt signal safely.
    """
    def __init__(self):
        super().__init__()
        self.emitter = QLogSignalEmitter()
        self.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))

    def emit(self, record):
        msg = self.format(record)
        self.emitter.log_signal.emit(msg)

# Global instance to attach to root/core loggers
qt_log_handler = QLogHandler()
logging.getLogger().addHandler(qt_log_handler)
logging.getLogger().setLevel(logging.INFO)
