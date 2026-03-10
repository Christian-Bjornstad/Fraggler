"""
Fraggler Diagnostics — PyQt6 Architecture Stylesheet
"""

VIBRANT_PRO_QSS = """
/* Global Application Settings */
QWidget {
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    font-size: 13px;
    color: #0f172a;
}

QMainWindow {
    background-color: #f8fafc;
}

/* Sidebar */
#Sidebar {
    background-color: #0f172a;
    border-right: 1px solid #1e293b;
}

#SidebarBrand {
    color: #ffffff;
    font-size: 16px;
    font-weight: 800;
    padding: 24px 20px 20px 20px;
    letter-spacing: 0.5px;
}

#SidebarButton {
    background: transparent;
    color: #94a3b8;
    text-align: left;
    padding: 12px 20px;
    border: none;
    font-weight: 600;
    font-size: 14px;
}

#SidebarButton:hover {
    background: #1e293b;
    color: #ffffff;
}

#SidebarButton:checked {
    background: #1e293b;
    color: #38bdf8;
    border-left: 4px solid #38bdf8;
    font-weight: 700;
}

/* Analysis Group Styles */
#AnalysisGroupHeader {
    background: transparent;
    color: #f8fafc;
    text-align: left;
    padding: 12px 20px;
    border: none;
    font-weight: 700;
    font-size: 14px;
    letter-spacing: 0.5px;
}

#AnalysisGroupHeader:hover {
    background: #1e293b;
}

#AnalysisGroupHeader:checked {
    color: #38bdf8;
    background: #1e293b;
}

#AnalysisSubButton {
    background: transparent;
    color: #94a3b8;
    text-align: left;
    padding: 10px 20px 10px 40px;
    border: none;
    font-weight: 500;
    font-size: 13px;
}

#AnalysisSubButton:hover {
    color: #ffffff;
    background: #1e293b;
}

#AnalysisSubButton:checked {
    color: #38bdf8;
    font-weight: 700;
}

/* Cards */
#Card {
    background: #ffffff;
    border-radius: 12px;
    border: 1px solid #e2e8f0;
}

#CardTitle {
    color: #64748b;
    font-size: 12px;
    font-weight: 800;
    text-transform: uppercase;
    letter-spacing: 1.2px;
    padding: 16px 16px 4px 16px;
}

/* Headers */
#PageTitle {
    font-size: 24px;
    font-weight: 800;
    color: #0f172a;
    line-height: 1.2;
}

#PageSubtitle {
    color: #64748b;
    font-size: 14px;
    margin-bottom: 20px;
}

/* Standard Buttons */
QPushButton {
    background-color: #ffffff;
    border: 1px solid #cbd5e1;
    border-radius: 8px;
    padding: 8px 18px;
    color: #334155;
    font-weight: 600;
    font-size: 13px;
}

QPushButton:hover {
    background-color: #f8fafc;
    border-color: #94a3b8;
    color: #0f172a;
}

QPushButton:pressed {
    background-color: #f1f5f9;
}

/* Primary Button */
QPushButton#PrimaryButton {
    background-color: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #06b6d4, stop:1 #4f46e5);
    border: none;
    color: #ffffff;
    border-radius: 8px;
    padding: 10px 24px;
    font-weight: 700;
    letter-spacing: 0.5px;
}

QPushButton#PrimaryButton:hover {
    background-color: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #0891b2, stop:1 #4338ca);
}

QPushButton#PrimaryButton:pressed {
    background-color: #3730a3;
}

QPushButton#PrimaryButton:disabled {
    background-color: #e2e8f0;
    color: #94a3b8;
}

/* Inputs */
QLineEdit, QDoubleSpinBox, QSpinBox {
    background-color: #ffffff;
    border: 1px solid #cbd5e1;
    border-radius: 8px;
    padding: 8px 12px;
    color: #0f172a;
    font-weight: 500;
    selection-background-color: #e0e7ff;
}

QLineEdit:focus, QDoubleSpinBox:focus, QSpinBox:focus {
    border: 1px solid #6366f1;
    background-color: #ffffff;
}

/* ComboBoxes (Premium Styling) */
QComboBox {
    background-color: #ffffff;
    border: 1px solid #cbd5e1;
    border-radius: 8px;
    padding: 8px 32px 8px 14px;
    color: #0f172a;
    font-weight: 600;
}

QComboBox:hover {
    border-color: #94a3b8;
}

QComboBox:focus {
    border: 1px solid #4f46e5;
}

QComboBox::drop-down {
    subcontrol-origin: padding;
    subcontrol-position: top right;
    width: 28px;
    border-left: none;
}

QComboBox QAbstractItemView {
    border: 1px solid #cbd5e1;
    border-radius: 8px;
    background-color: #ffffff;
    selection-background-color: #eff6ff;
    selection-color: #4f46e5;
    padding: 4px;
    outline: none;
}

/* Tabulator / Table styling */
QTableWidget {
    background-color: #ffffff;
    border: 1px solid #e2e8f0;
    border-radius: 8px;
    gridline-color: #f1f5f9;
}

QHeaderView::section {
    background-color: #f8fafc;
    color: #475569;
    font-weight: 700;
    text-transform: uppercase;
    border: none;
    border-bottom: 2px solid #e2e8f0;
    border-right: 1px solid #f1f5f9;
    padding: 6px;
    font-size: 11px;
}

QTableWidget::item:selected {
    background-color: #eff6ff;
    color: #0f172a;
}

/* Scrollbars */
QScrollBar:vertical {
    border: none;
    background: transparent;
    width: 10px;
    margin: 0;
}
QScrollBar::handle:vertical {
    background: #cbd5e1;
    min-height: 20px;
    border-radius: 5px;
}
QScrollBar::handle:vertical:hover {
    background: #94a3b8;
}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
    height: 0px;
}
"""
