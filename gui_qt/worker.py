import sys
import traceback
from PyQt6.QtCore import QObject, QRunnable, pyqtSignal, pyqtSlot

class WorkerSignals(QObject):
    '''
    Defines the signals available from a running worker thread.
    Supported signals are:
    
    finished: No data
    error: tuple (exctype, value, traceback.format_exc() )
    result: object data returned from processing
    progress: int indicating % progress
    status: str indicating status message
    log: str indicating log content to append
    '''
    finished = pyqtSignal()
    error = pyqtSignal(tuple)
    result = pyqtSignal(object)
    progress = pyqtSignal(int)
    progress_max = pyqtSignal(int)
    status = pyqtSignal(str)
    log = pyqtSignal(str)
    event = pyqtSignal(object)
    
    # Custom signal for batch jobs: (idx, total, name, state)
    progress_ext = pyqtSignal(int, int, str, str)

class Worker(QRunnable):
    '''
    Worker thread
    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.
    '''

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        # Add the callback to our kwargs if the target function requests it
        if 'progress_callback' in self.kwargs:
            self.kwargs['progress_callback'] = self.signals.progress
        if 'status_callback' in self.kwargs:
            self.kwargs['status_callback'] = self.signals.status
        if 'log_callback' in self.kwargs:
            self.kwargs['log_callback'] = self.signals.log
        if 'progress_max_callback' in self.kwargs:
            self.kwargs['progress_max_callback'] = self.signals.progress_max

    @pyqtSlot()
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''
        try:
            result = self.fn(*self.args, **self.kwargs)
        except Exception:
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)  # Return the result of the processing
        finally:
            self.signals.finished.emit()  # Done
