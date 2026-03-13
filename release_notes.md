# Fraggler Diagnostics v1.1.1 Release Notes

## New in v1.1.1
- **Fixed Basepair Shift**: Improved ladder fitting robustness by increasing the peak consideration limit and tightening contaminant filtering. This resolves observed shifts in samples with high primer/dye blobs (e.g., SL assay).

## v1.1.0 Overview
This version establishes Fraggler Diagnostics as a comprehensive, standalone solution for the automated processing and clinical validation of DNA fragment analysis data. The application is designed to provide clinical genomics labs with a secure, offline environment for high-throughput FSA trace analysis.

## Core Capabilities

Fraggler Diagnostics automates the traditionally manual and time-consuming steps of fragment analysis, from peak detection to quality control.

### Automated Clinical Pipeline
The tool features an intelligent scanning engine that automatically identifies assay types and categorizes datasets as either Quality Control (QC) or Patient samples. It performs peak alignment, fits size ladders, and applies passing/failing quality ratings based on validated laboratory metrics.

### High-Throughput Batch Processing
Labs can process multiple datasets simultaneously through the unified Run module. The system supports multi-folder input and drag-and-drop operations, building a centralized queue for background execution while providing real-time progress monitoring.

### Interactive Visual Diagnostics
The software generates high-resolution, interactive graphical documents (DIT Reports). These allows for fine-grained inspection of traces, allowing users to zoom into specific peaks, verify automated calls, and ensure data integrity within a diagnostic-grade interface.

### Longitudinal Quality Control
Fraggler maintains a centralized master trend log in Excel format. This enables labs to monitor ladder drift and assay performance over time, facilitating long-term quality assurance and longitudinal analysis.

### Clinical Reporting Standards
Output documentation is automatically formatted to meet clinical office standards. Reports are named according to established naming conventions and organized into date-based subdirectories for efficient archiving and retrieval.

## Security and Deployment

### Secure Offline Architecture
Fraggler Diagnostics is built for sensitive clinical environments. All data processing occurs locally on the host machine; no data is transmitted externally. The application bundles all necessary visualization assets and libraries, ensuring full functionality in air-gapped or restricted network environments.

### Cross-Platform Distribution
Standalone, pre-compiled binaries are now available for all major desktop platforms, requiring no local Python installation or complex setup:
- macOS (Native .app bundle)
- Windows (Standalone .exe)
- Linux (Portable ELF binary)

---
Developed for Oslo University Hospital (OUS) and the clinical genomics community.
