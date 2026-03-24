"""
Fraggler Diagnostics — PyQt6 Architecture Stylesheet
"""

VIBRANT_PRO_QSS = """
/* Global Application Settings */
QWidget {
    font-family: 'Avenir Next', 'Inter', 'Segoe UI', 'Helvetica Neue', Arial;
    font-size: 13px;
    color: #102235;
}

QMainWindow {
    background-color: #eef4f8;
}

QScrollArea#TabScrollArea {
    border: none;
    background: transparent;
}

QScrollArea#TabScrollArea > QWidget > QWidget {
    background: transparent;
}

/* Sidebar */
#Sidebar {
    background-color: #0b1724;
    border-right: 1px solid #17324a;
}

#SidebarBrand {
    color: #f8fbff;
    font-size: 17px;
    font-weight: 800;
    padding: 26px 20px 18px 20px;
    letter-spacing: 1.4px;
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

#SidebarButton:focus {
    background: #102235;
    color: #ffffff;
    border-left: 3px solid #7dd3fc;
    padding-left: 17px;
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
    color: #eff6ff;
    text-align: left;
    padding: 12px 20px;
    border: none;
    font-weight: 700;
    font-size: 14px;
    letter-spacing: 0.8px;
}

#AnalysisGroupHeader:hover {
    background: #102235;
}

#AnalysisGroupHeader:focus {
    background: #102235;
    border: 1px solid #7dd3fc;
}

#AnalysisGroupHeader:checked {
    color: #7dd3fc;
    background: #102235;
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
    background: #102235;
}

#AnalysisSubButton:focus {
    color: #ffffff;
    background: #122a3f;
    border-left: 3px solid #93c5fd;
    padding-left: 37px;
}

#AnalysisSubButton:checked {
    color: #7dd3fc;
    font-weight: 700;
    background: #122a3f;
    border-left: 3px solid #60a5fa;
    padding-left: 37px;
}

/* Cards */
#Card {
    background: #fbfdff;
    border-radius: 18px;
    border: 1px solid #d8e5ef;
}

#CardTitle {
    color: #486177;
    font-size: 12px;
    font-weight: 800;
    text-transform: uppercase;
    letter-spacing: 1.2px;
    padding: 18px 18px 6px 18px;
}

#DashboardCard {
    background: qlineargradient(
        x1: 0, y1: 0, x2: 1, y2: 1,
        stop: 0 #fdfefe,
        stop: 0.55 #f3f8fc,
        stop: 1 #e8f1f7
    );
    border: 1px solid #cfe0eb;
    border-radius: 20px;
}

#DashboardTitle {
    color: #0f2539;
    font-size: 23px;
    font-weight: 850;
}

#DashboardMetricCard {
    background: #ffffff;
    border: 1px solid #d6e2ec;
    border-radius: 14px;
}

#DashboardMetricLabel {
    color: #587185;
    font-size: 11px;
    font-weight: 800;
    text-transform: uppercase;
    letter-spacing: 1.2px;
}

#DashboardMetricValue {
    color: #0f2539;
    font-size: 19px;
    font-weight: 800;
}

#DashboardMetricDetail {
    color: #597285;
    font-size: 11px;
    line-height: 1.2;
}

#WorkflowSummaryText {
    color: #537083;
    font-size: 12px;
    font-weight: 700;
}

#WorkflowStatusBadge {
    border-radius: 14px;
    padding: 10px 14px;
    font-size: 12px;
    font-weight: 900;
    letter-spacing: 1.1px;
}

#WorkflowStatusBadge[state="ready"] {
    color: #1e4b71;
    background: #eaf4fb;
    border: 1px solid #bed8ec;
}

#WorkflowStatusBadge[state="running"] {
    color: #1d4f7a;
    background: #dcecff;
    border: 1px solid #9dc3ea;
}

#WorkflowStatusBadge[state="success"] {
    color: #16613c;
    background: #dcf4e6;
    border: 1px solid #9fd0b3;
}

#WorkflowStatusBadge[state="warning"] {
    color: #8a4f08;
    background: #fff3d8;
    border: 1px solid #ebcb8b;
}

#WorkflowStatusBadge[state="error"] {
    color: #8c2c1f;
    background: #ffe2de;
    border: 1px solid #ebb2a9;
}

#WorkflowStatusText {
    font-size: 14px;
    font-weight: 700;
}

#WorkflowStatusText[state="ready"] {
    color: #496375;
}

#WorkflowStatusText[state="running"] {
    color: #205b89;
}

#WorkflowStatusText[state="success"] {
    color: #1f7a4d;
}

#WorkflowStatusText[state="warning"] {
    color: #99580d;
}

#WorkflowStatusText[state="error"] {
    color: #b53b2d;
}

/* Headers */
#PageTitle {
    font-size: 27px;
    font-weight: 800;
    color: #0f2539;
    line-height: 1.2;
}

#PageSubtitle {
    color: #566f82;
    font-size: 14px;
    margin-bottom: 20px;
}

QLabel#MutedText {
    color: #6a8091;
}

/* Standard Buttons */
QPushButton {
    background-color: #ffffff;
    border: 1px solid #c9d8e4;
    border-radius: 10px;
    padding: 9px 18px;
    color: #23415a;
    font-weight: 600;
    font-size: 13px;
}

QPushButton:hover {
    background-color: #f2f7fb;
    border-color: #8fb7d4;
    color: #0f2539;
}

QPushButton:focus {
    border: 1px solid #2b6cb0;
    background-color: #f5faff;
}

QPushButton:pressed {
    background-color: #e9f1f8;
}

/* Primary Button */
QPushButton#PrimaryButton {
    background-color: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #2b6cb0, stop:1 #1d4f91);
    border: none;
    color: #ffffff;
    border-radius: 10px;
    padding: 10px 24px;
    font-weight: 700;
    letter-spacing: 0.5px;
}

QPushButton#PrimaryButton:hover {
    background-color: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #255f9b, stop:1 #18457f);
}

QPushButton#PrimaryButton:focus {
    border: 1px solid #0f3f72;
    background-color: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #255f9b, stop:1 #18457f);
}

QPushButton#PrimaryButton:pressed {
    background-color: #173b69;
}

QPushButton#PrimaryButton:disabled {
    background-color: #d9e3eb;
    color: #88a0b1;
}

/* Inputs */
QLineEdit, QDoubleSpinBox, QSpinBox {
    background-color: #ffffff;
    border: 1px solid #c9d8e4;
    border-radius: 10px;
    padding: 8px 12px;
    color: #0f2539;
    font-weight: 500;
    selection-background-color: #dcecff;
}

QLineEdit:focus, QDoubleSpinBox:focus, QSpinBox:focus {
    border: 1px solid #2b6cb0;
    background-color: #ffffff;
}

/* ComboBoxes (Premium Styling) */
QComboBox {
    background-color: #ffffff;
    border: 1px solid #c9d8e4;
    border-radius: 10px;
    padding: 8px 32px 8px 14px;
    color: #0f2539;
    font-weight: 600;
}

QComboBox:hover {
    border-color: #8fb7d4;
}

QComboBox:focus {
    border: 1px solid #2b6cb0;
}

QComboBox::drop-down {
    subcontrol-origin: padding;
    subcontrol-position: top right;
    width: 28px;
    border-left: none;
}

QComboBox QAbstractItemView {
    border: 1px solid #c9d8e4;
    border-radius: 10px;
    background-color: #ffffff;
    selection-background-color: #e7f1fb;
    selection-color: #1d4f91;
    padding: 4px;
    outline: none;
}

QFrame#GeneralSelectorCard,
QFrame#GeneralTraceCard {
    background: #f5f9fc;
    border: 1px solid #d6e3ed;
    border-radius: 14px;
}

QFrame#GeneralTraceCard:hover {
    border-color: #8eb7d3;
    background: #f7fbff;
}

QFrame#GeneralTraceCard[checked="true"] {
    border-color: #2b6cb0;
    background: #eaf3fb;
}

QCheckBox#GeneralTraceCheckbox {
    font-size: 13px;
    font-weight: 700;
    color: #0f2539;
    spacing: 10px;
}

QCheckBox#GeneralTraceCheckbox::indicator {
    width: 18px;
    height: 18px;
}

/* Tabulator / Table styling */
QTableWidget {
    background-color: #ffffff;
    border: 1px solid #d8e4ec;
    border-radius: 14px;
    gridline-color: #edf4f8;
    alternate-background-color: #f7fbfd;
}

QHeaderView::section {
    background-color: #eff5f8;
    color: #456176;
    font-weight: 700;
    text-transform: uppercase;
    border: none;
    border-bottom: 1px solid #d5e2eb;
    border-right: 1px solid #e6eef4;
    padding: 8px 6px;
    font-size: 11px;
}

QTableWidget::item:selected {
    background-color: #e7f1fb;
    color: #0f2539;
}

QListWidget {
    background: #ffffff;
    border: 1px solid #d8e4ec;
    border-radius: 14px;
    padding: 6px;
}

QListWidget::item {
    padding: 8px 10px;
    border-radius: 8px;
}

QListWidget::item:selected {
    background: #e7f1fb;
    color: #0f2539;
}

QProgressBar {
    min-height: 12px;
    max-height: 12px;
    border-radius: 6px;
    background: #dce8f0;
    border: none;
    text-align: center;
    color: transparent;
}

QProgressBar::chunk {
    border-radius: 6px;
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #2b6cb0, stop:1 #4ea6d8);
}

/* Scrollbars */
QScrollBar:vertical {
    border: none;
    background: transparent;
    width: 10px;
    margin: 0;
}
QScrollBar::handle:vertical {
    background: #b2c5d3;
    min-height: 20px;
    border-radius: 5px;
}
QScrollBar::handle:vertical:hover {
    background: #7f99ae;
}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
    height: 0px;
}
"""
