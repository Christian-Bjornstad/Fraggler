from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QLabel, QPlainTextEdit, QPushButton, QHBoxLayout
)
from PyQt6.QtCore import pyqtSignal, pyqtSlot
from PyQt6.QtGui import QFont

class TabLog(QWidget):
    # We can emit to this signal from anywhere if we get a reference to the Global TabLog,
    # or just use it within the UI thread directly
    log_signal = pyqtSignal(str)
    
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        
        # Header
        header_layout = QHBoxLayout()
        title = QLabel("System Log")
        title.setObjectName("PageTitle")
        sub = QLabel("Real-time execution logs from the core Fraggler pipeline.")
        sub.setObjectName("PageSubtitle")
        
        self.clear_btn = QPushButton("Clear Log")
        self.clear_btn.clicked.connect(self.clear_log)
        
        header_vbox = QVBoxLayout()
        header_vbox.addWidget(title)
        header_vbox.addWidget(sub)
        
        header_layout.addLayout(header_vbox)
        header_layout.addStretch()
        header_layout.addWidget(self.clear_btn)
        
        # Text Editor
        self.text_edit = QPlainTextEdit()
        self.text_edit.setReadOnly(True)
        # Monospaced font for logs
        font = QFont("Menlo", 12)
        font.setStyleHint(QFont.StyleHint.Monospace)
        self.text_edit.setFont(font)
        self.text_edit.setStyleSheet("background-color: #1e293b; color: #e2e8f0; padding: 8px; border-radius: 6px;")
        
        layout.addLayout(header_layout)
        layout.addWidget(self.text_edit)
        
        self.log_signal.connect(self.append_log)
        
    @pyqtSlot(str)
    def append_log(self, text: str):
        self.text_edit.appendPlainText(text)
        # Auto-scroll to bottom
        scrollbar = self.text_edit.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())
        
    def clear_log(self):
        self.text_edit.clear()
