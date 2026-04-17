from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QLabel, QTextBrowser, QVBoxLayout, QWidget

from gui_qt.about_content import (
    APP_OVERVIEW,
    THIRD_PARTY_NOTICE_PATH,
    THIRD_PARTY_SOFTWARE,
    UPSTREAM_LICENSE_PATH,
    load_text,
)


class TabAbout(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(18)

        header = QVBoxLayout()
        title = QLabel(APP_OVERVIEW["title"])
        title.setObjectName("PageTitle")
        subtitle = QLabel(APP_OVERVIEW["subtitle"])
        subtitle.setObjectName("PageSubtitle")
        subtitle.setWordWrap(True)
        header.addWidget(title)
        header.addWidget(subtitle)
        layout.addLayout(header)

        layout.addWidget(self._build_summary_card())
        layout.addWidget(self._build_third_party_card())
        layout.addWidget(self._build_notice_card())
        layout.addWidget(self._build_license_card())
        layout.addStretch()

    def _build_summary_card(self) -> QWidget:
        card = QWidget()
        card.setObjectName("Card")
        body = QVBoxLayout(card)
        body.setContentsMargins(18, 18, 18, 18)
        body.setSpacing(10)

        version = QLabel(APP_OVERVIEW["version_label"])
        version.setStyleSheet("font-weight: 700; color: #0f2539;")
        owner = QLabel(APP_OVERVIEW["owner_context"])
        owner.setWordWrap(True)
        status = QLabel(APP_OVERVIEW["repo_license_status"])
        status.setObjectName("MutedText")
        status.setWordWrap(True)

        body.addWidget(version)
        body.addWidget(owner)
        body.addWidget(status)
        return card

    def _build_third_party_card(self) -> QWidget:
        card = QWidget()
        card.setObjectName("Card")
        body = QVBoxLayout(card)
        body.setContentsMargins(18, 18, 18, 18)
        body.setSpacing(8)

        title = QLabel("Third-party software")
        title.setObjectName("CardTitle")
        title.setContentsMargins(0, 0, 0, 0)
        body.addWidget(title)

        text = QTextBrowser()
        text.setOpenExternalLinks(True)
        text.setMinimumHeight(170)
        text.setReadOnly(True)
        text.setMarkdown(self._third_party_markdown())
        body.addWidget(text)
        return card

    def _build_notice_card(self) -> QWidget:
        card = QWidget()
        card.setObjectName("Card")
        body = QVBoxLayout(card)
        body.setContentsMargins(18, 18, 18, 18)
        body.setSpacing(8)

        title = QLabel("Repository notice")
        title.setObjectName("CardTitle")
        title.setContentsMargins(0, 0, 0, 0)
        body.addWidget(title)

        note = QTextBrowser()
        note.setOpenExternalLinks(True)
        note.setReadOnly(True)
        note.setMinimumHeight(210)
        note.setMarkdown(
            f"Full notice file: `{THIRD_PARTY_NOTICE_PATH.name}`\n\n{load_text(THIRD_PARTY_NOTICE_PATH)}"
        )
        body.addWidget(note)
        return card

    def _build_license_card(self) -> QWidget:
        card = QWidget()
        card.setObjectName("Card")
        body = QVBoxLayout(card)
        body.setContentsMargins(18, 18, 18, 18)
        body.setSpacing(8)

        title = QLabel("Upstream MIT license")
        title.setObjectName("CardTitle")
        title.setContentsMargins(0, 0, 0, 0)
        body.addWidget(title)

        path_label = QLabel(f"License file: {UPSTREAM_LICENSE_PATH.name}")
        path_label.setAlignment(Qt.AlignmentFlag.AlignLeft)
        body.addWidget(path_label)

        license_text = QTextBrowser()
        license_text.setReadOnly(True)
        license_text.setMinimumHeight(220)
        license_text.setPlainText(load_text(UPSTREAM_LICENSE_PATH))
        body.addWidget(license_text)
        return card

    def _third_party_markdown(self) -> str:
        chunks: list[str] = []
        for item in THIRD_PARTY_SOFTWARE:
            chunks.append(
                "\n".join(
                    [
                        f"**{item['name']}**",
                        f"- Homepage: {item['homepage']}",
                        f"- Authors: {item['authors']}",
                        f"- Copyright: {item['copyright']}",
                        f"- License: {item['license_name']}",
                        f"- Summary: {item['summary']}",
                        f"- Local derived paths: {', '.join(item['derived_paths'])}",
                    ]
                )
            )
        return "\n\n".join(chunks)
