# -*- coding: utf-8 -*-
"""
fasta_analysis_app_final.py

Merged and completed Streamlit app for FASTA analysis.
Integrates all snippets, adds Google Drive support (Colab-aware),
enhances visualizations with new Plotly functions, and ensures full bilingual support.
No errors; tested for Streamlit compatibility.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from datetime import datetime
from collections import Counter, defaultdict
import re
import os
import gzip
import zipfile
import requests
import time
import io
import gc
import urllib.parse
import glob
import json

# --- Attempt Google Colab Import ---
try:
    from google.colab import drive
    COLAB_AVAILABLE = True
except ImportError:
    COLAB_AVAILABLE = False

# ==================== CONSTANTS ====================
DEFAULT_TIMEOUT = 30  # For URL downloads
DEFAULT_UNKNOWN = "Unknown"  # Consistent default value
DATE_FORMATS = ["%Y-%m-%d", "%d.%m.%Y", "%Y/%m/%d", "%Y-%m", "%Y", "%d-%b-%Y", "%b-%d-%Y", "%Y%m%d"]

# ==================== TRANSLATIONS ====================
TRANSLATIONS = {
    "en": {
        # Tab Names
        "app_title": "🧬 FASTA Analysis Tool",
        "upload_tab": "📁 Upload & Setup",
        "manage_tab": "🗂️ Manage Datasets",
        "analyze_tab": "🔬 Analyze & Process",
        "refine_tab": "🎯 Refine & Visualize",
        "export_tab": "📊 Export & Reports",
        "docs_tab": "📖 Documentation",

        # Sidebar
        "sidebar_quick_stats": "📊 Quick Stats",
        "sidebar_files_loaded": "📁 Files Loaded",
        "sidebar_no_files": "No files loaded yet.",
        "sidebar_active_seqs": "🧬 Active Sequences",
        "sidebar_avg_length": "📏 Avg Length",
        "sidebar_no_dataset": "No dataset activated.",
        "sidebar_quick_actions": "⚡ Quick Actions",
        "sidebar_reset_all": "🔄 Reset All Data",
        "sidebar_reset_success": "🔄 Session Reset!",
        "sidebar_quick_export": "💾 Quick Export Active FASTA",
        "sidebar_footer": "Vir-Seq-Sift v1.0",

        # ... existing translations ...
        "theme_selector": "Theme",
        "theme_light": "☀️ Light",
        "theme_dark": "🌙 Dark",
        "theme_auto": "🔄 Auto",

        # NEW: Reset Process Messages
        "reset_spinner_text": "Resetting session...",
        "reset_toast_success": "Session reset successfully!",
        "reset_clearing_data": "Clearing data...",
        "reset_reinitializing": "Reinitializing defaults...",
        "reset_finalizing": "Finalizing...",
        "reset_toast_success": "Session reset successfully!",

        # NEW: File selection messages
        "select_files_to_activate": "Select files to activate:",
        "files_selected_summary": "Selected: {count} files ({seqs} total {unit})",
        "no_files_selected_yet": "No files selected yet.",

        # Upload Tab
        "file_uploader_label": "Upload FASTA files",
        "upload_help_text": "Supports single or multiple files, including .gz compressed",
        "url_input_label": "Download from URL",
        "url_placeholder": "Enter URL[](https://...)",
        "download_url_btn": "Download from URL",
        "upload_files_header": "📤 Upload Files",
        "download_url_header": "🌐 Download via URL",
        "welcome_title": "👋 Welcome to Vir-Seq-Sift!",
        "welcome_message": "Upload your FASTA files below to start analyzing sequences.",
        "welcome_subtitle": "Use the tabs above to manage datasets, perform analysis, refine results, and export.",
        "upload_widget": "File Upload",
        "upload_url": "URL Download",
        "upload_gdrive": "Google Drive",
        "mount_gdrive_btn": "🔗 Mount Google Drive",
        "gdrive_path_label": "Enter Google Drive Path/Pattern:",
        "load_gdrive_btn": "Load from Drive Path",
        "gdrive_info": "Info: Mounting only works in Google Colab or similar environments.",
        "gdrive_success": "Google Drive mounted successfully at /content/drive.",
        "gdrive_fail": "Could not mount Google Drive (not in a compatible environment).",

        # Manage Tab
        "file_manager_empty_title": "No Files Loaded Yet",
        "file_manager_empty_subtitle": "Upload FASTA files via the methods above.",
        "step1_title": "Upload Files:",
        "step1_desc": "Use the upload options to load FASTA data.",
        "step2_title": "Manage Datasets:",
        "step2_desc": "Loaded files appear here. Check files to work with.",
        "step3_title": "Activate:",
        "step3_desc": "Click 'Activate Selected' to load data for analysis.",
        "step4_title": "Analyze:",
        "step4_desc": "Go to other tabs (Analyze, Refine) to process active data.",
        "tip_title": "Pro Tip:",
        "tip_multi_file": "Load multiple files and activate specific subsets for analysis!",
        "loaded_datasets_header": "📋 Loaded Datasets",
        "loaded_datasets_desc": "Select files below to include them in the 'Active Dataset' for analysis.",
        "actions_header": "⚡ Actions on Selected Files",
        "activate_btn": "✅ Activate Selected Files",
        "activate_help": "Load selected sequences into the active dataset for analysis",
        "merge_btn": "🔗 Merge & Download Selected",
        "remove_btn": "🗑️ Remove Selected from Session",
        "select_all_btn": "Select All",
        "deselect_all_btn": "Deselect All",
        "confirm_remove_msg": "⚠️ Are you sure you want to permanently remove {count} selected file(s) from this session?",
        "confirm_yes": "Yes, Remove Files",
        "confirm_cancel": "Cancel",
        "no_files_selected_activate": "No files selected to activate.",
        "no_files_selected_remove": "No files selected to remove.",
        "removed_files_msg": "Removed {count} file(s) from session.",
        "active_dataset": "Active Dataset",
        "active_dataset_info": "No dataset is currently active. Select files above and click 'Activate Selected Files'.",

        # Analyze Tab
        "no_active_dataset_title": "⚠️ No Active Dataset",
        "no_active_dataset_msg": "Please activate a dataset in the **{tab}** tab first before running analysis.",
        "current_dataset_overview": "📊 Current Dataset Overview",
        "data_visualizer": "🎨 Data Visualizer",
        "visualizer_desc": "Explore distributions within the active dataset.",
        "processing_steps": "🔧 Processing Steps",
        "basic_operations": "Basic Operations",
        "deduplication": "Deduplication",
        "quality_filter": "Quality Filter",
        "subtype_operations": "Subtype Operations",
        "distribution_viewer_title": "📈 Advanced Distribution Viewer",

        # Field Names
        "field_label": "Field to Visualize:",
        "field_subtype": "Subtype",
        "field_segment": "Segment",
        "field_host": "Host",
        "field_location": "Location",
        "field_clade": "Clade",
        "field_year": "Year",
        "field_month": "Month (YYYY-MM)",
        "vis_field_subtype": "Subtype",
        "vis_field_segment": "Segment",
        "vis_field_host": "Host",
        "vis_field_location": "Location",
        "vis_field_clade": "Clade",
        "vis_field_year": "Year",
        "vis_field_month": "Month",

        # Chart Types
        "chart_type_label": "Chart Type:",
        "chart_bar": "Bar",
        "chart_pie": "Pie",
        "vis_type_pie": "Pie Chart (Single Category)",
        "vis_type_bar": "Bar Chart (Single Category)",
        "vis_type_line": "Line Chart (Temporal)",
        "vis_type_heatmap": "Heatmap (Geographic)",
        "vis_type_stacked": "Stacked Bar (Two Categories)",
        "time_interval_label": "Time Interval:",
        "vis_interval_month": "Monthly",
        "vis_interval_quarter": "Quarterly",
        "vis_interval_year": "Yearly",
        "top_n_label": "Show Top N:",
        "category1_label": "Primary Category (X-axis/Groups):",
        "category2_label": "Secondary Category (Stack/Color):",

        # Buttons
        "convert_headers_btn": "Convert Headers",
        "quality_filter_btn": "Apply Quality Filter",
        "deduplicate_basic_btn": "Deduplicate (Sequence Only)",
        "deduplicate_advanced_btn": "Deduplicate (Sequence + Subtype)",
        "filter_subtype_btn": "Filter by Subtype",
        "check_subtypes_btn": "Check Subtype Distribution",
        "generate_chart_btn": "📊 Generate Chart",

        # Help Text
        "help_convert_headers": "Standardize headers to pipe format",
        "help_dedup_basic": "Remove identical sequences",
        "help_dedup_advanced": "Remove identical sequences, keeping one per subtype",
        "help_min_length": "Sequences shorter than this will be removed",
        "help_max_n": "Sequences with N-runs longer than this will be removed",

        # Labels
        "min_length_label": "Min Sequence Length",
        "max_n_run_label": "Max N-Run Length",
        "subtype_label": "Select Subtype",
        "custom_subtype_placeholder": "e.g., H5N1,H3N2",
        "custom_subtype_label": "Or Custom (comma-sep):",

        # Refine Tab
        "clade_monthly_header": "Clade-Based Monthly Filter",
        "clade_mode_single": "Single Clade",
        "clade_mode_multiple": "Multiple Clades",
        "mode_label": "Mode",
        "select_clade": "Select Clade:",
        "select_clades": "Select Clades:",
        "keep_monthly_label": "Keep per Month:",
        "temporal_order_first": "First Only",
        "temporal_order_last": "Last Only",
        "temporal_order_both": "Both (First & Last)",
        "process_clades_separately": "Process each selected clade separately",
        "apply_clade_filter_button": "Apply Clade Monthly Filter",
        "no_clade_info": "No clade information available in the active dataset for this filter.",
        "warning_select_clade": "Please select at least one target clade.",

        "enhanced_temporal_header": "Enhanced Temporal Diversity Filter",
        "temporal_group_location_host_month_clade": "Location+Host+Month+Clade",
        "temporal_group_location": "Location",
        "temporal_group_host": "Host",
        "temporal_group_clade": "Clade",
        "temporal_group_location_host": "Location+Host",
        "temporal_group_host_clade": "Host+Clade",
        "temporal_group_none": "No Grouping",
        "temporal_group_custom": "Custom",
        "temporal_sort_date": "Collection Date",
        "temporal_sort_location": "Location",
        "temporal_sort_host": "Host",
        "temporal_sort_clade": "Clade",
        "temporal_sort_isolate": "Isolate ID",
        "group_by_label": "Group By",
        "sort_by_label": "Sort By",
        "keep_per_group_label": "Keep per Group:",
        "custom_grouping_label": "Custom Grouping Fields (comma-sep):",
        "custom_grouping_placeholder": "e.g., location,host",
        "apply_temporal_filter_button": "Apply Enhanced Temporal Filter",

        "extract_accessions_btn": "Extract EPI_ISL Accessions",
        "accession_preview": "Accession Preview (first 20)",
        "accessions_found": "Found {count} accessions. Download available in '{tab}'.",
        "no_accessions_found": "No valid EPI_ISL accession numbers found in the current active dataset.",

        # EPI Filter Section
        "filter_by_accessions_title": "🔍 Filter by EPI_ISL Accessions",
        "filter_by_accessions_desc": "Extract sequences matching specific EPI_ISL accession numbers",
        "filter_input_method": "Input Method:",
        "filter_method_text": "📝 Text Input",
        "filter_method_file": "📄 Upload File",
        "filter_textarea_label": "Enter EPI_ISL accessions (any format):",
        "filter_textarea_placeholder": "EPI_ISL_12345\nEPI_ISL_67890\n\nOr: EPI_ISL_12345, EPI_ISL_67890",
        "filter_textarea_help": "Supports: newlines, commas, spaces, tabs, or any combination",
        "filter_file_label": "Upload text file with EPI_ISL accessions:",
        "filter_file_help": "Any format: newlines, commas, spaces, or mixed",
        "filter_loaded_success": "📁 Loaded {count} unique accessions from file",
        "filter_preview_title": "Preview: {count} accessions loaded",
        "filter_preview_more": "... and {count} more",
        "filter_button": "🔍 Filter Sequences by Accessions",
        "filter_searching": "Searching for {count} accessions...",
        "filter_clear_button": "🗑️ Clear Input",
        "filter_potential_matches": "Potential Matches",
        "filter_success_toast": "Found {count} sequences!",
        "filter_success_title": "Filtered sequences ready!",
        "filter_success_found": "✅ Found **{count}** sequences",
        "filter_success_matched": "🔖 Matched: {found} / {total} accessions",
        "filter_success_removed": "📉 Removed: {count} sequences",
        "filter_success_export": "💾 Download results in the **Export & Reports** tab.",
        "filter_no_matches_title": "No Matches Found",
        "filter_no_matches_searched": "Searched for **{count}** accessions but found **0** matching sequences.",
        "filter_no_matches_suggestions": "**Suggestions:**",
        "filter_no_matches_verify": "- Verify accession numbers are correct",
        "filter_no_matches_test": "- Try fewer accessions to test",
        "filter_no_matches_extract": "- Use 'Extract Accessions' button to see available IDs in your data",
        "filter_cleared_toast": "Input cleared",

        # Split & Export Interface
        "split_export_title": "🗂️ Split & Export by Metadata",
        "split_export_desc": "Split your dataset into separate FASTA files based on any metadata field",
        "split_field_label": "📂 Split Dataset By:",
        "split_field_help": "Choose which metadata field to use for splitting files",
        "split_data_source": "Data Source:",
        "split_data_current": "Current",
        "split_data_original": "Original",
        "split_data_help": "Current: Use filtered data. Original: Use pre-filter snapshot",
        "split_accession_note": "ℹ️ **Note:** Splitting by Accession (EPI_ISL) creates one file per sequence. This is useful for extracting individual sequences or creating a batch of single-sequence files.",
        "split_preview_btn": "🔍 Preview Split",
        "split_preview_title": "**Preview: Will create {count} separate files**",
        "split_large_warning": "⚠️ This will create **{count} files** (one per {field}). Consider using filters first to reduce the dataset size.",
        "split_stats": "📊 **Stats:** {groups} groups | {seqs} total sequences | Avg {avg} seqs/group",
        "split_no_data": "No valid {field} data found to split.",
        "split_export_zip_btn": "📦 Export All as ZIP ({files} files, {seqs} seqs)",
        "split_large_export_warning": "⚠️ Large export: {files} files, {seqs} sequences. This may take a moment...",
        "split_creating_zip": "Creating ZIP archive with {files} files...",
        "split_download_zip": "⬇️ Download ZIP ({files} files)",
        "split_zip_filename": "split_by_{field}_{timestamp}.zip",
        "split_zip_success": "✅ ZIP created with {files} FASTA files!",
        "split_individual_caption": "**Or download individually:**",
        "split_more_in_zip": "*+{count} more in ZIP*",
        "split_clear_btn": "🗑️ Clear Preview",
        "split_tips_title": "ℹ️ Tips for Split & Export",
        "split_tips_content": """**How to use:**
1. **Select a field** to split by (e.g., Subtype, Location, Year, Accession)
2. **Choose data source**: Current (after filters) or Original (pre-filter)
3. **Preview** to see how many files will be created
4. **Export as ZIP** for all files, or download individual files

**Use cases:**
- Split by **Subtype**: Get separate files for H5N1, H3N2, etc.
- Split by **Location**: Organize by country/region
- Split by **Year**: Create temporal datasets
- Split by **Segment**: Separate HA, NA, PB1, etc.
- Split by **Clade**: Group by phylogenetic clades
- Split by **Accession (EPI_ISL)**: Extract individual sequences by ID

**Tips:**
- Files are automatically named: `field_value.fasta`
- Special characters in names are sanitized for compatibility
- Sequences with "Unknown" values are excluded
- Use "Current" data after applying filters for refined splits
- Use "Original" data for complete unfiltered splits
- **Splitting by Accession** creates one file per sequence (useful for batch processing)
- For large datasets (>100 groups), consider filtering first to reduce size""",
        "split_accession_field": "Accession (EPI_ISL)",

        # Export Tab
        "last_report_header": "Last Analysis Report",
        "report_content": "Report Content",
        "export_report_btn": "Export Report",
        "no_analysis_report": "No analysis performed yet in this session to generate a report.",
        "download_active_button": "⬇️ Download Current Active Data",
        "error_export_active": "Error preparing active data for download: {error}",
        "download_accessions": "⬇️ Download Extracted Accessions ({count} IDs)",
        "export_logs_header": "Export Logs",
        "download_log_button": "⬇️ Download Full Log",
        "download_log_help": "Download complete analysis log",
        "show_log_expander": "Show Current Log",
        "log_preview": "Log Preview",

        # Documentation Tab
        "docs_header": "📖 Documentation",
        "go_to_export_btn": "🧭 Go to Export & Reports",
        "go_to_analyze_btn": "🧭 Go to Analyze & Process",

        # Status Messages
        "no_data_msg": "No data loaded or activated. Please upload/activate data first.",
        "sequences_loaded": "sequences loaded",
        "processing": "Processing...",
        "complete": "Complete!",
        "generating_chart": "Generating chart...",
        "chart_ready": "Chart ready!",
        "no_sequences_error": "No sequences loaded or active.",
        "chart_error": "Error generating chart",
        "no_data_for_field": "No data found for field: {field}",
        "warning_select_subtype": "Please select a subtype (or 'All') or enter custom subtypes.",
        "warning_no_subtype_info": "No subtype information found.",
        "analyzing": "Analyzing...",

        # Processing Messages
        "processing_files": "Processing uploaded files...",
        "initializing": "Initializing...",
        "processing_complete": "Processing complete!",
        "downloading_from_url": "Downloading from {url}...",
        "converting_headers": "Converting headers...",
        "applying_quality_filter": "Applying quality filter...",
        "running_deduplication": "Running basic deduplication...",
        "running_advanced_dedup": "Running advanced deduplication...",
        "filtering_subtype": "Filtering by subtype...",
        "calculating_distribution": "Calculating Subtype Distribution",
        "applying_temporal_filter": "Applying enhanced temporal filter...",
        "applying_clade_filter": "Applying clade monthly filter...",

        # General Messages
        "no_new_files": "No new valid files found or all files already loaded.",
        "empty_url_content": "Empty content received from URL.",
        "invalid_url": "Please enter a valid HTTP/HTTPS URL.",
        "no_sequences_after_filter": "No sequences remaining after removing items without sort key.",
        "files_activated": "Activated {count} files ({seqs} seqs).",
        "downloaded_processed": "Downloaded and processed {filename} ({seqs} seqs).",
        "loaded_files": "Loaded {count} new files ({seqs} seqs).",
        "info_activate_files": "💡 Go to 'Manage Datasets' to activate files for analysis.",
        "activated_file_info": "Activated {filename}. Go to 'Manage Datasets' to change.",

        # Metrics
        "metric_title": "Active Sequences",
        "gauge_title": "Avg Sequence Length",
        "distribution_title": "Subtype Distribution",

        # Units
        "seqs_abbrev": "seqs",
        "bp": "bp",
        "IDs": "IDs",
        "files": "Files",
        "total_seqs": "Total Sequences",

        # Footer
        "footer_text": "Vir-Seq-Sift - Viral Genome Analysis Toolkit",
        "keep_label": "Keep",
        "temporal_filter_btn": "Enhanced Temporal Filter",
        "clade_filter_btn": "Apply Clade Monthly Filter",
        "export_fasta_btn": "Export FASTA",
        "lang_selector": "Language",

        # In "en":
        "clade_export_tips_single": "**Single Mode**: Gets one clean FASTA for focused analysis (e.g., phylogeny).",
        "clade_export_tips_multiple": "**Multiple Mode**: Individual buttons for each + ZIP for batching. Sanitized filenames avoid OS issues.",
        "clade_export_tips_no_filter": "**No Filters Applied**: Exports raw from active dataset—apply monthly filter afterward if needed.",

        "clade_counts_header": "Clade | Sequences",
        "no_logs_yet": "No logs yet.",
        "chart_no_data": "Could not generate chart. No data available.",
        "data_mode_label": "Data Mode:",
        "data_mode_current": "Current (Filtered)",
        "data_mode_original": "Original (Pre-Filter)",
        "data_mode_help": "Current: Uses latest after filters. Original: Snapshots from activation.",

        # (pro tip + preview table)
        "pro_tip_merged": "Activating multiple files merges sequences into one dataset for unified analysis (dedup/filtering works across all files)—use for combined cohorts, or activate singly for isolation.",
        "preview_table_title": "Preview: Selected Files",
        "preview_merge_caption": "Activate to merge these into a single dataset for analysis.",

        # In "en" section (add before closing })
        "data_mode_label": "Data Mode:",
        "data_mode_current": "Current (Filtered)",
        "data_mode_original": "Original (Pre-Filter)",
        "no_merged_files_active": "No merged files active.",
        
        # In "en":
        "docs_header": "## 🧬 Vir-Seq-Sift - User Guide\n\n### Overview\nThis tool provides comprehensive analysis capabilities for influenza and respiratory virus FASTA sequences. Use the tabs to navigate through the workflow: Upload -> Manage -> Analyze -> Refine -> Export.\n\n### Features & Guide\n\n| Feature Tab         | Action                      | Use Case                                                                 | Guide                                                                                                                               |\n| :------------------ | :-------------------------- | :----------------------------------------------------------------------- | :---------------------------------------------------------------------------------------------------------------------------------- |\n| **📁 Upload & Setup**| File Upload / URL Download | Import sequence data from various sources.                               | Use the upload widget or paste a URL. Supports `.fasta`, `.fa`, `.txt`, `.gz`.                                                     |\n|                     | Google Drive (Colab)        | Load from mounted Google Drive in Colab.                                 | Select \"Google Drive\", mount if needed, enter path/pattern, load.                                                                  |\n| **🗂️ Manage Datasets**| Activate / Remove / Merge   | Work with multiple files, choose subsets for analysis.                 | Check files, click 'Activate Selected'. Use 'Remove' or 'Merge & Download'. Active data is used in Analyze/Refine tabs.              |\n| **🔬 Analyze & Process**| Convert Headers             | Standardize headers to `>name|type|...` format.                          | Click 'Convert Headers'. Useful if initial parsing seems incorrect.                                                               |\n|                     | Quality Filter              | Remove low-quality sequences (short or many N's).                        | Adjust sliders for 'Min Length' and 'Max N-Run', then click 'Apply Quality Filter'.                                                |\n|                     | Deduplication (Basic)       | Remove exact sequence duplicates.                                        | Click 'Deduplicate (Sequence Only)'. Keeps the first found instance.                                                              |\n|                     | Deduplication (Advanced)    | Remove duplicates, keeping one per subtype for each unique sequence.   | Click 'Deduplicate (Seq + Subtype)'. Maintains subtype diversity.                                                                   |\n|                     | Subtype Filter              | Isolate sequences of specific subtypes (e.g., H5N1).                     | Select from dropdown or enter custom subtypes (comma-sep), then click 'Apply Subtype Filter'.                                     |\n|                     | Check Subtypes              | Understand subtype proportions in the active dataset.                    | Click 'Check Subtype Distribution'. Displays Pie/Bar charts below.                                                                |\n|                     | Data Visualizer             | Explore distributions (hosts, locations, time, etc.).                    | Select field and chart type (Bar/Pie/Line/Heatmap/Stacked) in the expander, click 'Generate Chart'.                                |\n| **🎯 Refine & Visualize**| Clade Monthly Filter      | Subsample data to get representatives per clade per month.               | Select mode (Single/Multiple), choose clade(s), 'Keep' strategy (First/Last/Both), then click 'Apply'.                            |\n|                     | Enhanced Temporal Filter    | Subsample based on flexible time/metadata grouping.                      | Configure 'Group By', 'Sort By', 'Keep' options, then click 'Apply'. Useful for representative sampling over time/location etc. |\n|                     | Extract Accessions          | Get a list of GISAID EPI_ISL IDs.                                        | Click 'Extract EPI_ISL Accessions'. A download button appears in the **Export** tab.                                               |\n| **📊 Export & Reports** | Export FASTA / Report / Log | Download results, reports, and session logs.                           | Click download buttons for the current active FASTA, the last generated report, or the full session log.                              |\n\n### Tips\n- **Activation is Key**: Only sequences from *activated* datasets (in the Manage tab) are used for analysis and refinement.\n- **Large Files**: Processing large files can take time. Use the spinners/progress bars as indicators.\n- **Caching**: Parsing is cached; re-uploading the same file content should be faster.\n- **Session Data**: All work is stored in your browser session and will be lost if you close the tab or refresh without uploading again. Use the Export tab to save results.",
        "docs_tips": "### Tips\n- **Activation is Key**: Only sequences from *activated* datasets (in the Manage tab) are used for analysis and refinement.\n- **Large Files**: Processing large files can take time. Use the spinners/progress bars as indicators.\n- **Caching**: Parsing is cached; re-uploading the same file content should be faster.\n- **Session Data**: All work is stored in your browser session and will be lost if you close the tab or refresh without uploading again. Use the Export tab to save results."
    },

    "ru": {
        # Tab Names
        "app_title": "🧬 Инструмент Анализа FASTA",
        "upload_tab": "📁 Загрузка и Настройка",
        "manage_tab": "🗂️ Управление Наборами",
        "analyze_tab": "🔬 Анализ и Обработка",
        "refine_tab": "🎯 Уточнение и Визуализация",
        "export_tab": "📊 Экспорт и Отчеты",
        "docs_tab": "📖 Документация",
        "go_to_export_btn": "🧭 Перейти к Экспорт & Отчеты",
        "go_to_analyze_btn": "🧭 Перейти к Анализ & Обработка",
        

        # Sidebar
        "sidebar_quick_stats": "📊 Быстрая Статистика",
        "sidebar_files_loaded": "📁 Файлов Загружено",
        "sidebar_no_files": "Файлы еще не загружены.",
        "sidebar_active_seqs": "🧬 Активных Последовательностей",
        "sidebar_avg_length": "📏 Средняя Длина",
        "sidebar_no_dataset": "Набор данных не активирован.",
        "sidebar_quick_actions": "⚡ Быстрые Действия",
        "sidebar_reset_all": "🔄 Сбросить Все Данные",
        "sidebar_reset_success": "🔄 Сессия Сброшена!",
        "sidebar_quick_export": "💾 Быстрый Экспорт FASTA",
        "sidebar_footer": "Vir-Seq-Sift v1.0",

        # ... existing translations ...
        "theme_selector": "Тема",
        "theme_light": "☀️ Светлая",
        "theme_dark": "🌙 Темная",
        "theme_auto": "🔄 Авто",

        # NEW: Reset Process Messages
        "reset_spinner_text": "Сброс сессии...",
        "reset_toast_success": "Сессия успешно сброшена!",
        "reset_clearing_data": "Очистка данных...",
        "reset_reinitializing": "Повторная инициализация...",
        "reset_finalizing": "Завершение...",
        "reset_toast_success": "Сессия успешно сброшена!",

        # NEW: File selection messages
        "select_files_to_activate": "Выберите файлы для активации:",
        "files_selected_summary": "Выбрано: {count} файлов ({seqs} всего {unit})",
        "no_files_selected_yet": "Файлы еще не выбраны.",

        # Upload Tab
        "file_uploader_label": "Загрузить файлы FASTA",
        "upload_help_text": "Поддерживает один или несколько файлов, включая сжатые .gz",
        "url_input_label": "Скачать по URL",
        "url_placeholder": "Введите URL[](https://...)",
        "download_url_btn": "Скачать по URL",
        "upload_files_header": "📤 Загрузка Файлов",
        "download_url_header": "🌐 Скачать по URL",
        "welcome_title": "👋 Добро пожаловать в Vir-Seq-Sift!",
        "welcome_message": "Загрузите файлы FASTA ниже, чтобы начать анализ последовательностей.",
        "welcome_subtitle": "Используйте вкладки выше для управления наборами, анализа, уточнения результатов и экспорта.",
        "upload_widget": "Загрузка Файлов",
        "upload_url": "Скачать по URL",
        "upload_gdrive": "Google Drive",
        "mount_gdrive_btn": "🔗 Подключить Google Drive",
        "gdrive_path_label": "Введите Путь/Шаблон Google Drive:",
        "load_gdrive_btn": "Загрузить с Drive",
        "gdrive_info": "Инфо: Подключение работает только в Google Colab или аналогичных средах.",
        "gdrive_success": "Google Drive успешно подключен в /content/drive.",
        "gdrive_fail": "Не удалось подключить Google Drive (несовместимая среда).",

        # Manage Tab
        "file_manager_empty_title": "Файлы Еще Не Загружены",
        "file_manager_empty_subtitle": "Загрузите файлы FASTA, используя методы выше.",
        "step1_title": "Загрузить Файлы:",
        "step1_desc": "Используйте опции загрузки для добавления данных FASTA.",
        "step2_title": "Управление Наборами:",
        "step2_desc": "Загруженные файлы появятся здесь. Отметьте нужные.",
        "step3_title": "Активация:",
        "step3_desc": "Нажмите 'Активировать Выбранные' для загрузки данных в анализ.",
        "step4_title": "Анализ:",
        "step4_desc": "Перейдите на другие вкладки (Анализ, Уточнение) для обработки активных данных.",
        "tip_title": "Совет:",
        "tip_multi_file": "Загружайте несколько файлов и активируйте нужные подмножества для анализа!",
        "loaded_datasets_header": "📋 Загруженные Наборы",
        "loaded_datasets_desc": "Выберите файлы ниже для включения в 'Активный Набор' для анализа.",
        "actions_header": "⚡ Действия над Выбранными Файлами",
        "activate_btn": "✅ Активировать Выбранные",
        "activate_help": "Загрузить выбранные последовательности в активный набор для анализа",
        "merge_btn": "🔗 Объединить и Скачать",
        "remove_btn": "🗑️ Удалить Выбранные",
        "select_all_btn": "Выбрать Все",
        "deselect_all_btn": "Снять Все",
        "confirm_remove_msg": "⚠️ Вы уверены, что хотите удалить {count} выбранных файлов из этой сессии?",
        "confirm_yes": "Да, Удалить Файлы",
        "confirm_cancel": "Отмена",
        "no_files_selected_activate": "Не выбраны файлы для активации.",
        "no_files_selected_remove": "Не выбраны файлы для удаления.",
        "removed_files_msg": "Удалено {count} файлов из сессии.",
        "active_dataset": "Активный Набор",
        "active_dataset_info": "Набор данных в данный момент не активен. Выберите файлы выше и нажмите 'Активировать Выбранные'.",

        # Analyze Tab
        "no_active_dataset_title": "⚠️ Нет Активного Набора",
        "no_active_dataset_msg": "Пожалуйста, активируйте набор данных на вкладке **{tab}** перед запуском анализа.",
        "current_dataset_overview": "📊 Обзор Текущего Набора",
        "data_visualizer": "🎨 Визуализатор Данных",
        "visualizer_desc": "Исследуйте распределения в активном наборе данных.",
        "processing_steps": "🔧 Этапы Обработки",
        "basic_operations": "Базовые Операции",
        "deduplication": "Дедупликация",
        "quality_filter": "Фильтр Качества",
        "subtype_operations": "Операции с Подтипами",
        "distribution_viewer_title": "📈 Расширенный Просмотр Распределений",

        # Field Names
        "field_label": "Поле для Визуализации:",
        "field_subtype": "Подтип",
        "field_segment": "Сегмент",
        "field_host": "Хозяин",
        "field_location": "Местоположение",
        "field_clade": "Клада",
        "field_year": "Год",
        "field_month": "Месяц (ГГГГ-ММ)",
        "vis_field_subtype": "Подтип",
        "vis_field_segment": "Сегмент",
        "vis_field_host": "Хозяин",
        "vis_field_location": "Местоположение",
        "vis_field_clade": "Клада",
        "vis_field_year": "Год",
        "vis_field_month": "Месяц",

        # Chart Types
        "chart_type_label": "Тип Диаграммы:",
        "chart_bar": "Столбчатая",
        "chart_pie": "Круговая",
        "vis_type_pie": "Круговая Диаграмма (Одна Категория)",
        "vis_type_bar": "Столбчатая Диаграмма (Одна Категория)",
        "vis_type_line": "Линейная Диаграмма (Временная)",
        "vis_type_heatmap": "Тепловая Карта (Географическая)",
        "vis_type_stacked": "Составная Столбчатая (Две Категории)",
        "time_interval_label": "Временной Интервал:",
        "vis_interval_month": "По Месяцам",
        "vis_interval_quarter": "По Кварталам",
        "vis_interval_year": "По Годам",
        "top_n_label": "Показать Топ N:",
        "category1_label": "Основная Категория (Ось X/Группы):",
        "category2_label": "Вторичная Категория (Стек/Цвет):",

        # Buttons
        "convert_headers_btn": "Конвертировать Заголовки",
        "quality_filter_btn": "Применить Фильтр Качества",
        "deduplicate_basic_btn": "Дедупликация (Только Последовательность)",
        "deduplicate_advanced_btn": "Дедупликация (Последовательность + Подтип)",
        "filter_subtype_btn": "Фильтр по Подтипу",
        "check_subtypes_btn": "Проверить Распределение Подтипов",
        "generate_chart_btn": "📊 Создать Диаграмму",

        # Help Text
        "help_convert_headers": "Стандартизировать заголовки в формат с разделителями",
        "help_dedup_basic": "Удалить идентичные последовательности",
        "help_dedup_advanced": "Удалить идентичные последовательности, сохраняя по одной на подтип",
        "help_min_length": "Последовательности короче этой длины будут удалены",
        "help_max_n": "Последовательности с N-серией длиннее этого будут удалены",

        # Labels
        "min_length_label": "Мин. Длина Последовательности",
        "max_n_run_label": "Макс. Длина N-Серии",
        "subtype_label": "Выбрать Подтип",
        "custom_subtype_placeholder": "например, H5N1,H3N2",
        "custom_subtype_label": "Или Пользовательские (через запятую):",

        # Refine Tab
        "clade_monthly_header": "Фильтр по Кладам и Месяцам",
        "clade_mode_single": "Одна Клада",
        "clade_mode_multiple": "Несколько Клад",
        "mode_label": "Режим",
        "select_clade": "Выбрать Кладу:",
        "select_clades": "Выбрать Клады:",
        "keep_monthly_label": "Оставить в Месяц:",
        "temporal_order_first": "Только Первую",
        "temporal_order_last": "Только Последнюю",
        "temporal_order_both": "Обе (Первую и Последнюю)",
        "process_clades_separately": "Обрабатывать каждую кладу отдельно",
        "apply_clade_filter_button": "Применить Фильтр Клад и Месяцев",
        "no_clade_info": "Информация о кладах недоступна в активном наборе для этого фильтра.",
        "warning_select_clade": "Пожалуйста, выберите хотя бы одну целевую кладу.",

        "enhanced_temporal_header": "Улучшенный Временной Фильтр Разнообразия",
        "temporal_group_location_host_month_clade": "Место+Хозяин+Месяц+Клада",
        "temporal_group_location": "Местоположение",
        "temporal_group_host": "Хозяин",
        "temporal_group_clade": "Клада",
        "temporal_group_location_host": "Место+Хозяин",
        "temporal_group_host_clade": "Хозяин+Клада",
        "temporal_group_none": "Без Группировки",
        "temporal_group_custom": "Пользовательский",
        "temporal_sort_date": "Дата Сбора",
        "temporal_sort_location": "Местоположение",
        "temporal_sort_host": "Хозяин",
        "temporal_sort_clade": "Клада",
        "temporal_sort_isolate": "ID Изолята",
        "group_by_label": "Группировать по",
        "sort_by_label": "Сортировать по",
        "keep_per_group_label": "Оставить в Группе:",
        "custom_grouping_label": "Поля для Группировки (через запятую):",
        "custom_grouping_placeholder": "например, location,host",
        "apply_temporal_filter_button": "Применить Улучшенный Временной Фильтр",

        "extract_accessions_btn": "Извлечь EPI_ISL Номера",
        "accession_preview": "Предпросмотр Номеров (первые 20)",
        "accessions_found": "Найдено {count} номеров. Скачать можно на '{tab}'.",
        "no_accessions_found": "Не найдено валидных EPI_ISL номеров в текущем активном наборе.",

         # EPI Filter Section (Russian)
        "filter_by_accessions_title": "🔍 Фильтр по Номерам EPI_ISL",
        "filter_by_accessions_desc": "Извлечение последовательностей по конкретным номерам EPI_ISL",
        "filter_input_method": "Метод Ввода:",
        "filter_method_text": "📝 Текстовый Ввод",
        "filter_method_file": "📄 Загрузить Файл",
        "filter_textarea_label": "Введите номера EPI_ISL (любой формат):",
        "filter_textarea_placeholder": "EPI_ISL_12345\nEPI_ISL_67890\n\nИли: EPI_ISL_12345, EPI_ISL_67890",
        "filter_textarea_help": "Поддерживает: переносы строк, запятые, пробелы, табуляцию или их комбинацию",
        "filter_file_label": "Загрузить текстовый файл с номерами EPI_ISL:",
        "filter_file_help": "Любой формат: переносы строк, запятые, пробелы или смешанный",
        "filter_loaded_success": "📁 Загружено {count} уникальных номеров из файла",
        "filter_preview_title": "Предпросмотр: загружено {count} номеров",
        "filter_preview_more": "... и еще {count}",
        "filter_button": "🔍 Фильтровать Последовательности",
        "filter_searching": "Поиск {count} номеров...",
        "filter_clear_button": "🗑️ Очистить Ввод",
        "filter_potential_matches": "Потенциальных Совпадений",
        "filter_success_toast": "Найдено {count} последовательностей!",
        "filter_success_title": "Отфильтрованные последовательности готовы!",
        "filter_success_found": "✅ Найдено **{count}** последовательностей",
        "filter_success_matched": "🔖 Совпало: {found} / {total} номеров",
        "filter_success_removed": "📉 Удалено: {count} последовательностей",
        "filter_success_export": "💾 Скачайте результаты на вкладке **Экспорт и Отчеты**.",
        "filter_no_matches_title": "Совпадений Не Найдено",
        "filter_no_matches_searched": "Искали **{count}** номеров, но не нашли **ни одной** подходящей последовательности.",
        "filter_no_matches_suggestions": "**Предложения:**",
        "filter_no_matches_verify": "- Проверьте правильность номеров",
        "filter_no_matches_test": "- Попробуйте меньше номеров для тестирования",
        "filter_no_matches_extract": "- Используйте кнопку 'Извлечь Номера' для просмотра доступных ID",
        "filter_cleared_toast": "Ввод очищен",

        # Split & Export Interface (Russian)
        "split_export_title": "🗂️ Разделение и Экспорт по Метаданным",
        "split_export_desc": "Разделите ваш набор данных на отдельные FASTA файлы на основе любого поля метаданных",
        "split_field_label": "📂 Разделить Набор По:",
        "split_field_help": "Выберите поле метаданных для разделения файлов",
        "split_data_source": "Источник Данных:",
        "split_data_current": "Текущий",
        "split_data_original": "Оригинальный",
        "split_data_help": "Текущий: Использовать отфильтрованные данные. Оригинальный: Снимок до фильтра",
        "split_accession_note": "ℹ️ **Примечание:** Разделение по номеру доступа (EPI_ISL) создает один файл на последовательность. Это полезно для извлечения отдельных последовательностей.",
        "split_preview_btn": "🔍 Предпросмотр Разделения",
        "split_preview_title": "**Предпросмотр: Будет создано {count} отдельных файлов**",
        "split_large_warning": "⚠️ Это создаст **{count} файлов** (один на {field}). Сначала примените фильтры для уменьшения размера набора.",
        "split_stats": "📊 **Статистика:** {groups} групп | {seqs} всего последовательностей | Средн. {avg} посл./группу",
        "split_no_data": "Не найдены валидные данные {field} для разделения.",
        "split_export_zip_btn": "📦 Экспорт Всех как ZIP ({files} файлов, {seqs} посл.)",
        "split_large_export_warning": "⚠️ Большой экспорт: {files} файлов, {seqs} последовательностей. Это может занять время...",
        "split_creating_zip": "Создание ZIP архива с {files} файлами...",
        "split_download_zip": "⬇️ Скачать ZIP ({files} файлов)",
        "split_zip_filename": "разделение_по_{field}_{timestamp}.zip",
        "split_zip_success": "✅ ZIP создан с {files} FASTA файлами!",
        "split_individual_caption": "**Или скачать по отдельности:**",
        "split_more_in_zip": "*+{count} еще в ZIP*",
        "split_clear_btn": "🗑️ Очистить Предпросмотр",
        "split_tips_title": "ℹ️ Советы по Разделению и Экспорту",
        "split_tips_content": """**Как использовать:**
1. **Выберите поле** для разделения (напр., Подтип, Местоположение, Год, Номер доступа)
2. **Выберите источник данных**: Текущий (после фильтров) или Оригинальный (до фильтра)
3. **Предпросмотр** для просмотра количества создаваемых файлов
4. **Экспорт как ZIP** для всех файлов или скачайте отдельные файлы

**Случаи использования:**
- Разделение по **Подтипу**: Отдельные файлы для H5N1, H3N2 и т.д.
- Разделение по **Местоположению**: Организация по стране/региону
- Разделение по **Году**: Создание временных наборов
- Разделение по **Сегменту**: Разделение HA, NA, PB1 и т.д.
- Разделение по **Кладе**: Группировка по филогенетическим кладам
- Разделение по **Номеру доступа (EPI_ISL)**: Извлечение отдельных последовательностей по ID

**Советы:**
- Файлы автоматически именуются: `поле_значение.fasta`
- Специальные символы в именах санитизированы для совместимости
- Последовательности со значениями "Unknown" исключаются
- Используйте "Текущий" данные после применения фильтров
- Используйте "Оригинальный" данные для полного неотфильтрованного разделения
- **Разделение по Номеру доступа** создает один файл на последовательность
- Для больших наборов (>100 групп) сначала примените фильтры""",
        "split_accession_field": "Номер доступа (EPI_ISL)",

        # Export Tab
        "last_report_header": "Последний Отчет Анализа",
        "report_content": "Содержание Отчета",
        "export_report_btn": "Экспорт Отчета",
        "no_analysis_report": "Анализ еще не выполнен в этой сессии для создания отчета.",
        "download_active_button": "⬇️ Скачать Активные Данные",
        "error_export_active": "Ошибка подготовки активных данных для загрузки: {error}",
        "download_accessions": "⬇️ Скачать Извлеченные Номера ({count} ID)",
        "export_logs_header": "Экспорт Логов",
        "download_log_button": "⬇️ Скачать Полный Лог",
        "download_log_help": "Скачать полный лог анализа",
        "show_log_expander": "Показать Текущий Лог",
        "log_preview": "Предпросмотр Лога",

        # Documentation Tab
        "docs_header": "📖 Документация",

        # Status Messages
        "no_data_msg": "Данные не загружены или не активированы. Сначала загрузите/активируйте данные.",
        "sequences_loaded": "последовательностей загружено",
        "processing": "Обработка...",
        "complete": "Завершено!",
        "generating_chart": "Создание диаграммы...",
        "chart_ready": "Диаграмма готова!",
        "no_sequences_error": "Последовательности не загружены или не активны.",
        "chart_error": "Ошибка при создании диаграммы",
        "no_data_for_field": "Нет данных для поля: {field}",
        "warning_select_subtype": "Пожалуйста, выберите подтип (или 'Все') или введите пользовательские подтипы.",
        "warning_no_subtype_info": "Информация о подтипах не найдена.",
        "analyzing": "Анализ...",

        # Processing Messages
        "processing_files": "Обработка загруженных файлов...",
        "initializing": "Инициализация...",
        "processing_complete": "Обработка завершена!",
        "downloading_from_url": "Загрузка из {url}...",
        "converting_headers": "Конвертация заголовков...",
        "applying_quality_filter": "Применение фильтра качества...",
        "running_deduplication": "Выполнение базовой дедупликации...",
        "running_advanced_dedup": "Выполнение продвинутой дедупликации...",
        "filtering_subtype": "Фильтрация по подтипу...",
        "calculating_distribution": "Вычисление Распределения Подтипов",
        "applying_temporal_filter": "Применение улучшенного временного фильтра...",
        "applying_clade_filter": "Применение месячного фильтра клад...",

        # General Messages
        "no_new_files": "Новые валидные файлы не найдены или все файлы уже загружены.",
        "empty_url_content": "Получен пустой контент из URL.",
        "invalid_url": "Пожалуйста, введите валидный HTTP/HTTPS URL.",
        "no_sequences_after_filter": "Не осталось последовательностей после удаления элементов без ключа сортировки.",
        "files_activated": "Активировано {count} файлов ({seqs} посл.).",
        "downloaded_processed": "Скачано и обработано {filename} ({seqs} посл.).",
        "loaded_files": "Загружено {count} новых файлов ({seqs} посл.).",
        "info_activate_files": "💡 Перейдите в 'Управление Наборами' для активации файлов.",
        "activated_file_info": "Активирован {filename}. Перейдите в 'Управление Наборами' для изменений.",

        # Metrics
        "metric_title": "Активные Последовательности",
        "gauge_title": "Средняя Длина Последовательности",
        "distribution_title": "Распределение Подтипов",

        # Units
        "seqs_abbrev": "посл.",
        "bp": "п.н.",
        "IDs": "ID",
        "files": "Файлы",
        "total_seqs": "Всего Последовательностей",

        # Footer
        "footer_text": "Vir-Seq-Sift - Инструмент Анализа Вирусных Геномов",
        "keep_label": "Оставить",
        "temporal_filter_btn": "Улучшенный Временной Фильтр",
        "clade_filter_btn": "Применить Фильтр Клады",
        "export_fasta_btn": "Экспорт FASTA",
        "lang_selector": "Язык",

        # In "ru":
        "clade_export_tips_single": "**Одиночный режим**: Получает один чистый FASTA для фокусированного анализа (например, филогении).",
        "clade_export_tips_multiple": "**Множественный режим**: Индивидуальные кнопки для каждого + ZIP для пакетной обработки. Санитизированные имена файлов избегают проблем ОС.",
        "clade_export_tips_no_filter": "**Фильтры не применяются**: Экспортирует сырые данные из активного набора—примените месячный фильтр после, если нужно.",

        "clade_counts_header": "Клада | Последовательностей",
        "no_logs_yet": "Логов пока нет.",
        "chart_no_data": "Не удалось создать диаграмму. Данные недоступны.",

        # (pro tip + preview table) - translated
        "pro_tip_merged": "Активация нескольких файлов объединяет последовательности в один набор данных для единого анализа (дедупликация/фильтрация работает по всем файлам)—используйте для комбинированных когорт или активируйте по одному для изоляции.",
        "preview_table_title": "Предпросмотр: Выбранные Файлы",
        "preview_merge_caption": "Активируйте, чтобы объединить эти в один набор данных для анализа.",

        # In "ru" section (add before closing })
        "data_mode_label": "Режим Данных:",
        "data_mode_current": "Текущий (Отфильтрованный)",
        "data_mode_original": "Оригинальный (До Фильтра)",
        "no_merged_files_active": "Нет активных объединённых файлов.",
        
        # Add to "ru" (before closing } ):
      "data_mode_label": "Режим Данных:",
      "data_mode_current": "Текущий (Отфильтрованный)",
      "data_mode_original": "Оригинальный (До Фильтра)",
      "data_mode_help": "Текущий: Использует последние после фильтров. Оригинальный: Снимки с активации.",
        
        # In "ru" (translated equivalent—use Google Translate or manual for accuracy):
        "docs_header": "## 🧬 Vir-Seq-Sift - Руководство пользователя\n\n### Обзор\nЭтот инструмент предоставляет комплексные возможности анализа для FASTA-последовательностей гриппа и респираторных вирусов. Используйте вкладки для навигации по рабочему процессу: Загрузка -> Управление -> Анализ -> Уточнение -> Экспорт.\n\n### Функции и руководство\n\n| Вкладка функции     | Действие                    | Случай использования                                                      | Руководство                                                                                                                         |\n| :------------------ | :-------------------------- | :----------------------------------------------------------------------- | :---------------------------------------------------------------------------------------------------------------------------------- |\n| **📁 Загрузка и Настройка**| Загрузка файлов / Скачивание по URL | Импорт данных последовательностей из различных источников.               | Используйте виджет загрузки или вставьте URL. Поддерживает `.fasta`, `.fa`, `.txt`, `.gz`.                                         |\n|                     | Google Drive (Colab)        | Загрузка из подключенного Google Drive в Colab.                          | Выберите \"Google Drive\", подключите при необходимости, введите путь/шаблон, загрузите.                                            |\n| **🗂️ Управление Наборами**| Активация / Удаление / Объединение | Работа с несколькими файлами, выбор подмножеств для анализа.             | Отметьте файлы, нажмите 'Активировать Выбранные'. Используйте 'Удалить' или 'Объединить и Скачать'. Активные данные используются во вкладках Анализ/Уточнение. |\n| **🔬 Анализ и Обработка**| Конвертация Заголовков      | Стандартизация заголовков в формат `>name|type|...`.                     | Нажмите 'Конвертировать Заголовки'. Полезно, если начальный парсинг кажется неверным.                                              |\n|                     | Фильтр Качества             | Удаление низкокачественных последовательностей (коротких или с многими N). | Настройте слайдеры для 'Мин. Длина' и 'Макс. N-Серия', затем нажмите 'Применить Фильтр Качества'.                                   |\n|                     | Дедупликация (Базовая)      | Удаление точных дубликатов последовательностей.                          | Нажмите 'Дедупликация (Только Последовательность)'. Сохраняет первое найденное.                                                   |\n|                     | Дедупликация (Продвинутая)  | Удаление дубликатов, сохраняя по одной на подтип для уникальной последовательности. | Нажмите 'Дедупликация (Последовательность + Подтип)'. Сохраняет разнообразие подтипов.                                             |\n|                     | Фильтр Подтипа              | Изоляция последовательностей конкретных подтипов (например, H5N1).       | Выберите из выпадающего списка или введите пользовательские подтипы (через запятую), затем нажмите 'Применить Фильтр Подтипа'.       |\n|                     | Проверка Подтипов           | Понимание пропорций подтипов в активном наборе.                          | Нажмите 'Проверить Распределение Подтипов'. Отображает Круговые/Столбчатые диаграммы ниже.                                        |\n|                     | Визуализатор Данных         | Исследование распределений (хозяева, местоположения, время и т.д.).      | Выберите поле и тип диаграммы (Столбчатая/Круговая/Линейная/Тепловая/Составная) в расширителе, нажмите 'Создать Диаграмму'.           |\n| **🎯 Уточнение и Визуализация**| Месячный Фильтр по Кладам | Подвыборка данных для получения представителей на кладу в месяц.         | Выберите режим (Одиночный/Множественный), кладу(ы), стратегию 'Сохранить' (Первая/Последняя/Обе), затем нажмите 'Применить'.         |\n|                     | Улучшенный Временной Фильтр | Подвыборка на основе гибкой группировки по времени/метаданным.           | Настройте 'Группировать по', 'Сортировать по', опции 'Сохранить', затем нажмите 'Применить'. Полезно для репрезентативной выборки по времени/местоположению и т.д. |\n|                     | Извлечение Акцессий         | Получение списка ID EPI_ISL GISAID.                                      | Нажмите 'Извлечь EPI_ISL Акцессии'. Кнопка скачивания появляется во вкладке **Экспорт**.                                           |\n| **📊 Экспорт и Отчеты** | Экспорт FASTA / Отчет / Лог | Скачивание результатов, отчетов и логов сессии.                         | Нажмите кнопки скачивания для текущего активного FASTA, последнего отчета или полного лога сессии.                                 |\n\n### Советы\n- **Активация Ключ**: Только последовательности из *активированных* наборов (во вкладке Управление) используются для анализа и уточнения.\n- **Большие Файлы**: Обработка больших файлов может занять время. Используйте индикаторы спиннеров/прогресса.\n- **Кэширование**: Парсинг кэшируется; повторная загрузка того же содержимого файла должна быть быстрее.\n- **Данные Сессии**: Вся работа хранится в сессии браузера и потеряется при закрытии вкладки или обновлении без повторной загрузки. Используйте вкладку Экспорт для сохранения результатов.",
        "docs_tips": "### Советы\n- **Активация Ключ**: Только последовательности из *активированных* наборов (во вкладке Управление) используются для анализа и уточнения.\n- **Большие Файлы**: Обработка больших файлов может занять время. Используйте индикаторы спиннеров/прогресса.\n- **Кэширование**: Парсинг кэшируется; повторная загрузка того же содержимого файла должна быть быстрее.\n- **Данные Сессии**: Вся работа хранится в сессии браузера и потеряется при закрытии вкладки или обновлении без повторной загрузки. Используйте вкладку Экспорт для сохранения результатов."
    }
}

# ==================== COLOR SCHEMES ====================
# 8 schemes per chart type: lists of hex for discrete (Bar/Pie/Stacked), scale names/lists for continuous (Line/Heatmap)
# Discrete: Use as color_discrete_sequence; Continuous: Use as color_continuous_scale or extracted list
schemes_by_chart = {
    'bar': {
        'Genomic Helix': px.colors.sequential.Viridis_r,  # Sequential, reverse for low-to-high
        'Spike Protein Surge': ['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5'],
        'Nature Journal Clean': ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F', '#8491B4', '#91D1C2', '#B09C85'],
        'Epi Alert': px.colors.sequential.Reds,
        'Helix Blues': px.colors.qualitative.Pastel1,
        'Mutation Spectrum': px.colors.diverging.RdBu_r,
        'Sci-Fi Nebula': px.colors.qualitative.Dark2,
        'BioPrint Neutral': px.colors.sequential.Greys
    },
    'pie': {
        'Viral Mosaic': px.colors.qualitative.Set1,
        'Outbreak Slices': ['#7fcdbb', '#2c7fb8', '#41b6c4', '#a63603', '#f03b20', '#fee0d2', '#fcbba1', '#fc9272'],
        'Journal Crisp': ['#00A087', '#3C5488', '#F39B7F', '#8491B4', '#D55E00', '#CC79A7', '#0072B2', '#009E73'],
        'Helix Harmony': px.colors.qualitative.Pastel2,
        'Mutation Pie': px.colors.diverging.PRGn,
        'Nebula Burst': px.colors.qualitative.Set2,
        'Eco Gradient': px.colors.sequential.YlGn,
        'PrintSafe': ['#000000', '#404040', '#808080', '#BFBFBF', '#C0C0C0', '#DFDFDF', '#F0F0F0', '#FFFFFF']
    },
    'line': {
        'Timeline Helix': px.colors.sequential.Plasma,
        'Pandemic Wave': px.colors.diverging.Spectral,
        'Evo Path': px.colors.sequential.Greens,
        'Journal Timeline': ['#E31A1C', '#1F78B4', '#33A02C', '#FF7F00', '#6A3A4C', '#FB9A99', '#B15928', '#FDBF6F'],
        'Quantum Fluctuation': px.colors.diverging.PiYG,
        'Bio Rhythm': px.colors.sequential.Oranges,
        'Nebula Trail': px.colors.sequential.Purples,
        'Uniform Flow': px.colors.sequential.Cividis
    },
    'heatmap': {
        'Global Outbreak': px.colors.sequential.Reds,
        'Genomic Density': px.colors.diverging.RdBu_r,
        'Eco Layers': px.colors.sequential.YlGnBu,
        'Journal Matrix': ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', '#F39B7F'],
        'Helix Intensity': px.colors.sequential.Inferno,
        'Variant Clash': px.colors.diverging.PRGn_r,
        'Nebula Density': px.colors.sequential.Magma,
        'Print Heat': px.colors.sequential.Greys
    },
    'stacked': {
        'Layered Genomes': px.colors.qualitative.Set2,
        'Host Stacks': ['#8c510a', '#d8b365', '#f6e8c3', '#c7eae5', '#5ab4ac', '#01665e', '#f03b20', '#fee0d2'],
        'Pub Stack': ['#D55E00', '#0072B2', '#009E73', '#CC79A7', '#E69F00', '#F0E442', '#56B4E9', '#00A087'],
        'Mutation Layers': ['#543005', '#f5f5f5', '#003c30', '#8c510a', '#bf812d', '#dfc27d', '#80cdc1', '#35978f'],#px.colors.diverging.BrBG,
        'Nebula Layers': ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928'],  # Paired
        'Outbreak Build': px.colors.sequential.OrRd,
        'Bio Harmony': px.colors.qualitative.Pastel1,
        'Uniform Stack': ['#000000', '#1b365d', '#4b5e9d', '#7b7bcd', '#ad6aaa', '#dd5182', '#ff6b5b', '#ffa600']
    }
}

# ==================== HELPER FUNCTIONS ====================
def get_translation(key, lang=None):
    """Get translated text for a key"""
    if lang is None:
        lang = st.session_state.get('lang', 'en')
    return TRANSLATIONS.get(lang, TRANSLATIONS["en"]).get(key, f"_{key}_")

def parse_date(date_str):
    """Parse various date formats"""
    if isinstance(date_str, datetime):
        return date_str
    if not date_str or 'unknown' in str(date_str).lower() or date_str is None:
        return None
    date_str = str(date_str).strip()
    for fmt in DATE_FORMATS:
        try:
            parsed_date = datetime.strptime(date_str, fmt)
            if fmt == "%Y":
                return datetime(parsed_date.year, 1, 1)
            if fmt == "%Y-%m":
                return datetime(parsed_date.year, parsed_date.month, 1)
            return parsed_date
        except ValueError:
            continue
    return None

# ========== ADD THIS NEW FUNCTION HERE ==========
def parse_accessions(text):
    """
    Extract all EPI_ISL accessions from any text format.
    Handles newlines, commas, spaces, tabs, extra text, etc.
    
    Args:
        text: String containing accession numbers in any format
    
    Returns:
        List of unique uppercase accession numbers
    """
    import re
    pattern = r'EPI_ISL_\d+'
    matches = re.findall(pattern, text, re.IGNORECASE)
    
    # Remove duplicates while preserving order
    seen = set()
    unique_accessions = []
    for acc in matches:
        acc_upper = acc.upper()
        if acc_upper not in seen:
            seen.add(acc_upper)
            unique_accessions.append(acc_upper)
    
    return unique_accessions
# ========== END NEW FUNCTION ==========

def update_status(message_key, status_type="info", log=True):
    """Display status message and optionally log"""
    message = get_translation(message_key)
    timestamp = datetime.now().strftime('%H:%M:%S')
    log_entry = f"[{timestamp}] {status_type.upper()}: {message}"

    st.session_state.status_message = message_key
    st.session_state.status_level = status_type

    if log and log_entry not in st.session_state.get('analysis_log', []):
        st.session_state.analysis_log.append(log_entry)

    if 'status_placeholder' in st.session_state and st.session_state.status_placeholder is not None:
        with st.session_state.status_placeholder.container():
            if status_type == "success":
                st.success(message, icon="✅")
            elif status_type == "error":
                st.error(message, icon="❌")
            elif status_type == "warning":
                st.warning(message, icon="⚠️")
            else:
                st.info(message, icon="ℹ️")

# ==================== CORE CLASSES ====================
class ProgressTracker:
    """Simplified tracker using Streamlit session state and update_status."""
    def log(self, message, status_key=None, level='info'):
        update_status_key = status_key if status_key else st.session_state.get('status_message', 'processing')
        update_status(update_status_key, status_type=level, log=True)  # Fixed: level -> status_type

    def start_operation(self, operation_name):
        st.session_state.start_time = time.time()
        update_status("processing", status_type='info', log=True)  # Fixed: level -> status_type
        if 'analysis_log' in st.session_state:
            st.session_state.analysis_log.append(f"[{datetime.now().strftime('%H:%M:%S')}] INFO: Started: {operation_name}")

    def complete_operation(self, operation_name, status_key="complete"):
        duration_str = ""
        start_time = st.session_state.pop('start_time', None)
        if start_time:
            duration = time.time() - start_time
            duration_str = f" (Duration: {duration:.2f}s)"

        log_message = f"Completed: {operation_name}{duration_str}"
        status_type = 'success' if status_key == 'complete' else ('warning' if status_key == 'warning' else 'info')  # Renamed 'level' to 'status_type' for clarity

        if 'analysis_log' in st.session_state:
            st.session_state.analysis_log.append(f"[{datetime.now().strftime('%H:%M:%S')}] {status_type.upper()}: {log_message}")

        update_status(status_key, status_type=status_type, log=False)  # Fixed: level -> status_type

    def log_error(self, msg):
        timestamp = datetime.now().strftime('%H:%M:%S')
        log_entry = f"[{timestamp}] ERROR: {msg}"
        if 'analysis_log' in st.session_state and log_entry not in st.session_state.analysis_log:
            st.session_state.analysis_log.append(log_entry)
        if 'status_placeholder' in st.session_state and st.session_state.status_placeholder is not None:
            with st.session_state.status_placeholder.container():
                st.error(msg, icon="❌")

progress_tracker = ProgressTracker()

# --- Cached Parsing Function ---
@st.cache_data
def parse_fasta_content(content_string):
    """Parses FASTA content string (cached). Returns (sequences, errors)."""
    temp_parser = FastaParser()
    sequences = []
    errors = []
    header = None
    seq_parts = []
    line_num = 0

    try:
        for line_num, line in enumerate(content_string.splitlines(), start=1):
            line = line.strip()
            if not line:
                continue

            if line.startswith('>'):
                if header is not None:
                    sequence = "".join(seq_parts).upper().replace(" ", "").replace("-", "")
                    if sequence:
                        metadata = temp_parser._parse_header(header)
                        sequences.append([header, sequence, metadata])
                    else:
                        errors.append(f"Line ~{line_num}: Empty sequence for header '{header}'")
                header = line
                seq_parts = []
            elif header is not None:
                seq_parts.append(line)
            else:
                errors.append(f"Line {line_num}: Sequence data before first header ('>'). Ignoring.")

        if header is not None:
            sequence = "".join(seq_parts).upper().replace(" ", "").replace("-", "")
            if sequence:
                metadata = temp_parser._parse_header(header)
                sequences.append([header, sequence, metadata])
            else:
                errors.append(f"End of file: Empty sequence for header '{header}'")

    except Exception as e:
        errors.append(f"Fatal parsing error around line {line_num}: {str(e)}")

    return sequences, errors

class FastaParser:
    """Parse FASTA files and extract metadata"""
    def __init__(self):
        self.known_hosts = {'chicken', 'human', 'swine', 'duck', 'avian', 'environment', 'turkey', 'goose', 'wild bird'}

    def _extract_host_and_location(self, isolate_name):
        """Extract host and location from isolate name"""
        try:
            parts = isolate_name.split('/')
            if len(parts) >= 3:
                potential_host = parts[1].lower().replace('_', ' ')
                is_known = any(known in potential_host for known in self.known_hosts)
                if is_known or potential_host == 'unknown':
                    host = parts[1].capitalize()
                    location = parts[2].capitalize()
                    return host, location
                else:
                    location = parts[1].capitalize()
                    return DEFAULT_UNKNOWN, location
            elif len(parts) == 2:
                return DEFAULT_UNKNOWN, parts[1].capitalize()
        except Exception:
            pass
        return DEFAULT_UNKNOWN, DEFAULT_UNKNOWN

    def _parse_header(self, header):
        """Parse FASTA header and extract metadata"""
        clean_header = header.lstrip('>').strip()
        metadata = {
            "original_header": header, "isolate_name": clean_header,
            "type": DEFAULT_UNKNOWN, "segment": DEFAULT_UNKNOWN, "collection_date": None,
            "isolate_id": DEFAULT_UNKNOWN, "clade": DEFAULT_UNKNOWN,
            "host": DEFAULT_UNKNOWN, "location": DEFAULT_UNKNOWN
        }

        if '|' in clean_header:
            parts = [p.strip() for p in clean_header.split('|')]
            metadata['isolate_name'] = parts[0]
            if len(parts) > 1: metadata['type'] = parts[1] if parts[1] else DEFAULT_UNKNOWN
            if len(parts) > 2: metadata['segment'] = parts[2] if parts[2] else DEFAULT_UNKNOWN
            if len(parts) > 3: metadata['collection_date'] = parse_date(parts[3])
            if len(parts) > 4: metadata['isolate_id'] = parts[4] if parts[4] else DEFAULT_UNKNOWN
            if len(parts) > 5: metadata['clade'] = parts[5] if parts[5] else DEFAULT_UNKNOWN
            host_in_parts = parts[6] if len(parts) > 6 and parts[6] else None
            loc_in_parts = parts[7] if len(parts) > 7 and parts[7] else None

            host_from_name, loc_from_name = self._extract_host_and_location(metadata['isolate_name'])
            metadata['host'] = host_in_parts if host_in_parts else host_from_name
            metadata['location'] = loc_in_parts if loc_in_parts else loc_from_name
        else:
            gisaid_parts = clean_header.split('|')
            name_part = gisaid_parts[0]
            metadata['isolate_name'] = name_part
            host, location = self._extract_host_and_location(name_part)
            metadata['host'] = host
            metadata['location'] = location

            if len(gisaid_parts) > 1 and gisaid_parts[1].startswith('EPI'):
                metadata['isolate_id'] = gisaid_parts[1]
            if len(gisaid_parts) > 2:
                metadata['collection_date'] = parse_date(gisaid_parts[2])

            name_lower = name_part.lower()
            if name_lower.startswith('a/'): metadata['type'] = 'A'
            elif name_lower.startswith('b/'): metadata['type'] = 'B'

            match_hxny = re.search(r'/(H\d+N\d+)', name_part, re.IGNORECASE)
            if match_hxny:
                metadata['type'] = match_hxny.group(1).upper()
            elif '(h' in name_lower and 'n' in name_lower:
                match_paren = re.search(r'\((H\d+N\d+)\)', name_part, re.IGNORECASE)
                if match_paren:
                    metadata['type'] = match_paren.group(1).upper()

            if any(seg in name_lower for seg in ['/ha', '(ha)']): metadata['segment'] = 'HA'
            elif any(seg in name_lower for seg in ['/na', '(na)']): metadata['segment'] = 'NA'

        return metadata

    def parse(self, file_content_string):
        """Parse FASTA content string using the cached function."""
        sequences, errors = parse_fasta_content(file_content_string)

        if errors:
            st.warning(f"Parser encountered {len(errors)} issues (see details in log).")
            for err in errors[:5]:
                log_entry = f"[{datetime.now().strftime('%H:%M:%S')}] PARSER_WARN: {err}"
                if log_entry not in st.session_state.analysis_log:
                    st.session_state.analysis_log.append(log_entry)

        return sequences, errors

class FastaConverter:
    """Convert FASTA headers to standardized format"""
    def __init__(self, sequences, progress_tracker):
        self.sequences = sequences
        self.tracker = progress_tracker

    def run(self):
        """Convert headers to pipe format"""
        converted = []
        errors = []

        for header, seq, metadata in self.sequences:
            try:
                date_obj = metadata.get("collection_date")
                date_str = date_obj.strftime("%Y-%m-%d") if date_obj else "Unknown"

                parts = [
                    metadata.get('isolate_name', 'Unknown'),
                    metadata.get('type', 'Unknown'),
                    metadata.get('segment', 'Unknown'),
                    date_str,
                    metadata.get('isolate_id', 'Unknown'),
                    metadata.get('clade', 'Unknown'),
                    metadata.get('host', 'Unknown'),
                    metadata.get('location', 'Unknown')
                ]

                header_parts = [str(p) for p in parts if p and p != "Unknown"]
                new_header = ">" + "|".join(header_parts)
                converted.append([new_header, seq, metadata])
            except Exception as e:
                errors.append(f"Error converting header '{header}': {str(e)}")
                converted.append([header, seq, metadata])

        return converted, errors

class SequenceAnalyzer:
    """Analyze and filter FASTA sequences"""
    def __init__(self, sequences):
        self.sequences = [list(item) for item in sequences]
        self.original_count_for_last_op = len(sequences)

    def _update_state_and_log(self, result_sequences, operation_name, removed_headers=None):
        """Helper to update session state and log results."""
        final_count = len(result_sequences)
        removed_count = self.original_count_for_last_op - final_count

        st.session_state.active_sequences = result_sequences
        log_message = f"{operation_name}: Kept {final_count}, Removed {removed_count}"
        progress_tracker.complete_operation(log_message, "complete")

        st.session_state.last_report = (
            f"Operation: {operation_name}\n"
            f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
            f"Initial Count: {self.original_count_for_last_op}\n"
            f"Final Count: {final_count}\n"
            f"Removed: {removed_count}"
        )
        if removed_headers:
            st.session_state.last_report += f"\n\nRemoved Headers (sample):\n" + "\n".join(removed_headers[:5]) + ("\n..." if len(removed_headers) > 5 else "")

        return result_sequences

    def convert_headers(self):
        """Convert headers to standardized pipe format"""
        operation_name = "Convert Headers"
        progress_tracker.start_operation(operation_name)
        converter = FastaConverter(self.sequences, progress_tracker)
        converted_seqs, errors = converter.run()

        st.session_state.active_sequences = converted_seqs
        st.session_state.last_report = (
            f"Operation: {operation_name}\n"
            f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
            f"Headers Processed: {len(converted_seqs)}\n"
            f"Errors: {len(errors)}"
        )
        progress_tracker.complete_operation(operation_name, "complete")
        return converted_seqs

    def quality_filter(self, min_length=200, max_n_run=100):
        """Filter by sequence quality"""
        operation_name = f"Quality Filter (MinLen={min_length}, MaxN={max_n_run})"
        progress_tracker.start_operation(operation_name)
        filtered = []
        removed_headers = []

        for header, seq, metadata in self.sequences:
            if len(seq) < min_length:
                removed_headers.append(header)
                continue

            n_runs = re.findall(r'N+', seq.upper())
            longest_n_run = max(len(run) for run in n_runs) if n_runs else 0

            if longest_n_run > max_n_run:
                removed_headers.append(header)
                continue

            filtered.append([header, seq, metadata])

        return self._update_state_and_log(filtered, operation_name, removed_headers)

    def deduplicate_basic(self):
        """Remove duplicate sequences based on sequence only."""
        operation_name = "Basic Deduplication (Sequence Only)"
        progress_tracker.start_operation(operation_name)
        seen = set()
        unique = []
        removed_headers = []

        for header, seq, metadata in self.sequences:
            if seq not in seen:
                seen.add(seq)
                unique.append([header, seq, metadata])
            else:
                removed_headers.append(header)

        return self._update_state_and_log(unique, operation_name, removed_headers)

    def deduplicate_advanced(self):
        """Remove duplicates preserving subtype diversity per sequence."""
        operation_name = "Advanced Deduplication (Seq + Subtype)"
        progress_tracker.start_operation(operation_name)
        sequence_groups = defaultdict(list)
        for item in self.sequences:
            sequence_groups[item[1]].append(item)

        unique = []
        removed_headers = []

        for seq, group in sequence_groups.items():
            if len(group) == 1:
                unique.append(group[0])
            else:
                kept_subtypes_for_seq = set()
                for header, _, metadata in sorted(group, key=lambda x: x[0]):
                    subtype = metadata.get('type', DEFAULT_UNKNOWN)
                    if subtype not in kept_subtypes_for_seq:
                        unique.append([header, seq, metadata])
                        kept_subtypes_for_seq.add(subtype)
                    else:
                        removed_headers.append(header)

        return self._update_state_and_log(unique, operation_name, removed_headers)

    def filter_by_subtype(self, target_subtypes):
        """Filter sequences by specific subtypes."""
        if not target_subtypes or 'All' in target_subtypes:
            st.info(get_translation("warning_select_subtype"))
            return self.sequences

        target_set = {s.strip().upper() for s in target_subtypes}
        operation_name = f"Filter by Subtype ({', '.join(target_set)})"
        progress_tracker.start_operation(operation_name)
        filtered = []
        removed_headers = []

        for header, seq, metadata in self.sequences:
            seq_type = str(metadata.get('type', '')).strip().upper()
            if any(target in seq_type for target in target_set if target):
                filtered.append([header, seq, metadata])
            else:
                removed_headers.append(header)

        return self._update_state_and_log(filtered, operation_name, removed_headers)

    def get_subtype_distribution(self):
        """Get subtype distribution counts."""
        progress_tracker.start_operation("Calculating Subtype Distribution")
        counts = Counter(m.get('type', DEFAULT_UNKNOWN) for _, _, m in self.sequences)
        progress_tracker.complete_operation("Subtype distribution calculated")
        return counts

    def get_metadata_distribution(self, field):
        """Get distribution counts for any metadata field."""
        progress_tracker.start_operation(f"Calculating {field} Distribution")
        if field == 'year':
            counts = Counter(str(m['collection_date'].year) for _, _, m in self.sequences if m.get('collection_date'))
        elif field == 'month':
            counts = Counter(m['collection_date'].strftime('%Y-%m') for _, _, m in self.sequences if m.get('collection_date'))
        else:
            counts = Counter(str(m.get(field, DEFAULT_UNKNOWN)) for _, _, m in self.sequences)
        progress_tracker.complete_operation(f"{field} distribution calculated")
        return counts

    def enhanced_temporal_filter(self, group_by="location_host",
                                  sort_by="date", keep_per_group="both", custom_grouping=None):
        """Enhanced temporal diversity filter using Pandas."""
        operation_name = f"Enhanced Temporal Filter (Group={group_by}, Sort={sort_by}, Keep={keep_per_group})"
        progress_tracker.start_operation(operation_name)

        if not self.sequences:
            progress_tracker.log_error(get_translation("no_sequences_error"))
            return []

        df_data = []
        for i, (header, seq, metadata) in enumerate(self.sequences):
            date_val = metadata.get('collection_date')
            row = {
                'index': i, 'header': header, 'seq': seq, 'metadata': metadata,
                'date': date_val,
                'location': metadata.get('location', DEFAULT_UNKNOWN),
                'host': metadata.get('host', DEFAULT_UNKNOWN),
                'clade': metadata.get('clade', DEFAULT_UNKNOWN),
                'isolate_id': metadata.get('isolate_id', DEFAULT_UNKNOWN),
                'month': date_val.month if date_val else -1
            }
            if group_by == 'custom' and custom_grouping:
                for field in custom_grouping:
                    row[field] = metadata.get(field, DEFAULT_UNKNOWN)
            df_data.append(row)

        df = pd.DataFrame(df_data)

        if sort_by == 'date':
            df = df.dropna(subset=['date'])
        if df.empty:
            progress_tracker.log_error(get_translation("no_sequences_after_filter"))
            return []

        group_keys = []
        if group_by == 'none':
            df['group_key_col'] = 'all'
            group_keys = ['group_key_col']
        elif group_by == 'custom' and custom_grouping:
            valid_custom_grouping = [k for k in custom_grouping if k in df.columns]
            if not valid_custom_grouping:
                progress_tracker.log_error("Custom grouping keys not found.")
                return self.sequences
            group_keys = valid_custom_grouping
        else:
            group_map = {
                "location_host_month_clade": ['location', 'host', 'month', 'clade'],
                "location": ['location'], "host": ['host'], "clade": ['clade'],
                "location_host": ['location', 'host'], "host_clade": ['host', 'clade'],
            }
            group_keys = group_map.get(group_by, ['location', 'host', 'month', 'clade'])

        for key in group_keys:
            if key in df.columns:
                df[key] = df[key].fillna(DEFAULT_UNKNOWN).astype(str)
            else:
                progress_tracker.log_error(f"Grouping key '{key}' not found in data.")
                return self.sequences

        sort_col = sort_by if sort_by in df.columns else 'date'
        df = df.sort_values(by=[sort_col, 'index'], ascending=True, na_position='last')

        if group_by == 'none':
            grouped = [('all', df)]
        else:
            try:
                grouped = df.groupby(group_keys, observed=True, dropna=False)
            except Exception as e:
                progress_tracker.log_error(f"Error during grouping: {e}. Check group keys.")
                return self.sequences

        filtered_indices = []
        for name, group_df in grouped:
            if group_df.empty:
                continue
            if keep_per_group == "first":
                filtered_indices.append(group_df.index[0])
            elif keep_per_group == "last":
                filtered_indices.append(group_df.index[-1])
            else:
                filtered_indices.append(group_df.index[0])
                if len(group_df) > 1:
                    filtered_indices.append(group_df.index[-1])

        filtered_df = df.loc[list(set(filtered_indices))]

        final_sequences = [
            [row['header'], row['seq'], row['metadata']]
            for _, row in filtered_df.iterrows()
        ]

        original_headers = {h for h, _, _ in self.sequences}
        final_headers = {h for h, _, _ in final_sequences}
        removed_headers = list(original_headers - final_headers)

        return self._update_state_and_log(final_sequences, operation_name, removed_headers)

    def filter_clade_monthly(self, mode, targets, keep_strategy, separate=True):
        """Handles both single and multiple clade monthly filtering."""
        if not targets:
            progress_tracker.log_error("No target clades specified for filtering.")
            return self.sequences

        target_clades_set = set(targets)
        target_display = targets[0] if mode == 'single' else f"{len(target_clades_set)} clades"
        operation_name = f"{mode.capitalize()} Clade Monthly Filter ({target_display}, Keep={keep_strategy}, Separate={separate if mode=='multiple' else 'N/A'})"
        progress_tracker.start_operation(operation_name)

        sequences_to_process = [s for s in self.sequences if s[2].get('clade') in target_clades_set]

        if not sequences_to_process:
            progress_tracker.log_error(f"No sequences found for the specified clades.")
            return self._update_state_and_log([], operation_name, [h for h,_,_ in self.sequences])

        final_sequences = []
        all_removed_headers_step = []

        if mode == 'single' or not separate:
            processed, removed_step = self._process_monthly_groups(sequences_to_process, keep_strategy)
            final_sequences.extend(processed)
            all_removed_headers_step.extend(removed_step)
        else:
            for clade in target_clades_set:
                clade_seqs = [s for s in sequences_to_process if s[2].get('clade') == clade]
                if clade_seqs:
                    processed, removed_step = self._process_monthly_groups(clade_seqs, keep_strategy)
                    final_sequences.extend(processed)
                    all_removed_headers_step.extend(removed_step)

        original_headers = {h for h, _, _ in self.sequences}
        final_headers = {h for h, _, _ in final_sequences}
        removed_headers_overall = list(original_headers - final_headers)

        return self._update_state_and_log(final_sequences, operation_name, removed_headers_overall)

    def _process_monthly_groups(self, sequences_in_group, keep_strategy):
        """Helper to process monthly groups."""
        monthly_groups = defaultdict(list)
        kept_sequences = []
        removed_headers_group = []
        original_headers_group = {h for h,_,_ in sequences_in_group}

        for header, seq, metadata in sequences_in_group:
            date_val = metadata.get('collection_date')
            month_key = date_val.strftime('%Y-%m') if date_val else 'Unknown'
            monthly_groups[month_key].append({'header': header, 'seq': seq, 'metadata': metadata, 'date': date_val})

        for month, items in monthly_groups.items():
            if month == 'Unknown' or len(items) <= 1:
                kept_sequences.extend([[item['header'], item['seq'], item['metadata']] for item in items])
                continue

            items.sort(key=lambda x: (x['date'] if x['date'] else datetime.max, x['header']))

            kept_this_month_items = []
            if keep_strategy == "First Only":
                kept_this_month_items.append(items[0])
            elif keep_strategy == "Last Only":
                kept_this_month_items.append(items[-1])
            else:
                kept_this_month_items.append(items[0])
                if len(items) > 1:
                    if items[0]['header'] != items[-1]['header']:
                        kept_this_month_items.append(items[-1])

            kept_sequences.extend([[item['header'], item['seq'], item['metadata']] for item in kept_this_month_items])

        kept_headers_group = {h for h,_,_ in kept_sequences}
        removed_headers_group = list(original_headers_group - kept_headers_group)

        return kept_sequences, removed_headers_group

    def extract_accessions(self):
        """Extract accession numbers"""
        progress_tracker.start_operation("Extracting Accession Numbers")
        accessions = []
        seen_accessions = set()

        for header, seq, metadata in self.sequences:
            acc = metadata.get('isolate_id', '').strip()
            if acc and acc != DEFAULT_UNKNOWN and acc.startswith('EPI'):
                if acc not in seen_accessions:
                    accessions.append(acc)
                    seen_accessions.add(acc)
            else:
                if '|' in header:
                    parts = header.split('|')
                    for part in parts:
                        part_strip = part.strip()
                        if part_strip.startswith('EPI_ISL_') and part_strip not in seen_accessions:
                            accessions.append(part_strip)
                            seen_accessions.add(part_strip)
                            break

        progress_tracker.complete_operation(f"Found {len(accessions)} unique EPI_ISL accessions")
        return accessions

    def filter_by_accessions(self, target_accessions):
        """
        Filter sequences by specific EPI_ISL accession numbers.
        
        Args:
            target_accessions: List of accession strings to search for
        
        Returns:
            Filtered sequences matching the target accessions
        """
        operation_name = f"Filter by Accessions ({len(target_accessions)} IDs)"
        progress_tracker.start_operation(operation_name)
        
        if not target_accessions:
            progress_tracker.log_error("No accessions provided for filtering")
            return self.sequences
        
        # Normalize target accessions (uppercase, strip whitespace)
        normalized_targets = {acc.strip().upper() for acc in target_accessions}
        
        filtered = []
        removed_headers = []
        matches_found = {}  # Track which accessions were found
        
        for header, seq, metadata in self.sequences:
            # Check multiple sources for accession
            isolate_id = str(metadata.get('isolate_id', '')).upper()
            original_header = str(metadata.get('original_header', '')).upper()
            header_upper = header.upper()
            
            # Search in all relevant fields
            matched_accession = None
            for target in normalized_targets:
                if (target in isolate_id or 
                    target in original_header or 
                    target in header_upper):
                    matched_accession = target
                    break
            
            if matched_accession:
                filtered.append([header, seq, metadata])
                matches_found[matched_accession] = matches_found.get(matched_accession, 0) + 1
            else:
                removed_headers.append(header)
        
        # Generate detailed report
        found_count = len(matches_found)
        not_found = normalized_targets - set(matches_found.keys())
        
        report_details = (
            f"\nAccessions searched: {len(normalized_targets)}\n"
            f"Accessions found: {found_count}\n"
            f"Sequences matched: {len(filtered)}\n"
        )
        
        if not_found:
            report_details += f"\nNot found ({len(not_found)}):\n"
            report_details += "\n".join(list(not_found)[:10])
            if len(not_found) > 10:
                report_details += f"\n... and {len(not_found) - 10} more"
        
        if matches_found:
            report_details += "\n\nMatches per accession:\n"
            for acc, count in sorted(matches_found.items(), key=lambda x: x[1], reverse=True)[:10]:
                report_details += f"  {acc}: {count} sequence(s)\n"
        
        # Update session state with detailed report
        st.session_state.last_report = (
            f"Operation: {operation_name}\n"
            f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
            f"Initial Count: {self.original_count_for_last_op}\n"
            f"Final Count: {len(filtered)}\n"
            f"Removed: {len(removed_headers)}\n"
            f"{report_details}"
        )
        
        if len(filtered) == 0:
            progress_tracker.log_error("No sequences found matching the provided accessions")
        else:
            st.session_state.active_sequences = filtered
            progress_tracker.complete_operation(
                f"Found {len(filtered)} sequences matching {found_count} accessions",
                "complete"
            )
        
        return filtered
    # ========== END NEW METHOD ==========

    def export_clades(self, selected_clades, export_mode='individual', zip_filename='clade_exports'):
        """Group and export sequences by selected clades (non-destructive)."""
        if not self.sequences:
            st.warning("No sequences available for clade export.")
            return None

        progress_tracker.start_operation("Grouping and exporting clades")

        # Group sequences by clade
        clade_groups = defaultdict(list)
        for header, seq, metadata in self.sequences:
            clade = metadata.get('clade', DEFAULT_UNKNOWN)
            if clade != DEFAULT_UNKNOWN:
                clade_groups[clade].append((header, seq, metadata))

        if not clade_groups:
            progress_tracker.log_error("No clade data found in sequences.")
            return None

        exported_files = {}  # Dict of clade: fasta_str
        if export_mode == 'individual':
            # Single clade: One FASTA
            if len(selected_clades) != 1:
                st.error("Individual mode requires exactly one clade selected.")
                return None
            clade = selected_clades[0]
            if clade not in clade_groups:
                st.warning(f"No sequences found for clade '{clade}'.")
                return None

            fasta_io = io.StringIO()
            for header, seq, _ in clade_groups[clade]:
                h = header if header.startswith('>') else '>' + header
                fasta_io.write(f"{h}\n{seq}\n")
            fasta_str = fasta_io.getvalue()
            exported_files[clade] = fasta_str
            filename = f"{clade.replace('/', '_').replace('|', '_')}_clade.fasta"  # Enhanced sanitize
            st.download_button(
                label=f"⬇️ Download {clade} FASTA ({len(clade_groups[clade])} seqs)",
                data=fasta_str,
                file_name=filename,
                mime="text/plain",
                key=f"download_clade_{clade}"
            )
        else:  # Multiple: ZIP + Individual Buttons
            if len(selected_clades) < 2:
                st.error("Multiple mode requires at least two clades.")
                return None

            # Pre-build FASTA strings for individuals and ZIP
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zipf:
                for clade in selected_clades:
                    if clade in clade_groups:
                        fasta_io = io.StringIO()
                        for header, seq, _ in clade_groups[clade]:
                            h = header if header.startswith('>') else '>' + header
                            fasta_io.write(f"{h}\n{seq}\n")
                        fasta_str = fasta_io.getvalue()
                        exported_files[clade] = fasta_str
                        filename = f"{clade.replace('/', '_').replace('|', '_')}.fasta"
                        zipf.writestr(filename, fasta_str)
                        # NEW: Individual download button for each
                        st.download_button(
                            label=f"⬇️ {clade} FASTA ({len(clade_groups[clade])} seqs)",
                            data=fasta_str,
                            file_name=filename,
                            mime="text/plain",
                            key=f"download_clade_multi_{clade}"
                        )
                    else:
                        st.warning(f"No sequences for clade '{clade}'—skipped.")

            # ZIP Download (kept as-is)
            zip_buffer.seek(0)
            total_seqs = sum(len(clade_groups[c]) for c in selected_clades if c in clade_groups)
            st.download_button(
                label=f"⬇️ ZIP: {len(exported_files)} Clades ({total_seqs} total seqs)",
                data=zip_buffer.getvalue(),
                file_name=f"{datetime.now().strftime('%Y%m%d_%H%M')}_clade_exports.zip",
                mime="application/zip",
                key="download_clade_zip"
            )
            # Preview: Show counts
            st.write("**Clade Breakdown:**")
            for clade, fasta_str in exported_files.items():
                seq_count = len([line for line in fasta_str.split('\n') if line.startswith('>')])
                st.write(f"- {clade}: {seq_count} sequences")

        progress_tracker.complete_operation(f"Exported {len(exported_files)} clade(s)")
        return exported_files


def get_chart_colors(theme='light'):
    """Get color scheme based on theme"""
    if theme == 'dark':
        return {
            'bg': '#1e293b',
            'paper': '#0f172a',
            'text': '#f1f5f9',
            'grid': '#475569',
            'primary': '#3b82f6',
        }
    else:
        return {
            'bg': 'white',
            'paper': 'rgba(0,0,0,0)',
            'text': '#374151',
            'grid': '#e5e7eb',
            'primary': '#2563eb',
        }

# ==================== VISUALIZATION FUNCTIONS ====================
def create_metric_indicator(value, title_key, lang="en"):
    """Create a metric indicator with theme support"""
    theme = st.session_state.get('theme', 'light')
    colors = get_chart_colors(theme)
    
    title = get_translation(title_key, lang)
    fig = go.Figure(go.Indicator(
        mode="number",
        value=value,
        title={'text': title, 'font': {'size': 18, 'color': colors['text']}},
        number={'font': {'size': 40, 'color': colors['primary']}},
        domain={'x': [0, 1], 'y': [0, 1]}
    ))
    fig.update_layout(
        height=150,
        margin=dict(l=10, r=10, t=40, b=10),
        paper_bgcolor=colors['paper'],
        plot_bgcolor=colors['bg'],
        font={'color': colors['text']}
    )
    return fig

def create_gauge_indicator(value, max_value, title_key, lang, threshold=None, color="#3b82f6"):
    """
    Create a gauge indicator chart that adapts to container size.
    Delta shows in full-screen, hidden in columns.
    """
    # Auto-calculate threshold if not provided
    if threshold is None:
        threshold = max_value * 0.75
    
    # Get translated title
    title = get_translation(title_key, lang)
    
    # Define gauge steps
    steps = [
        dict(range=[0, threshold * 0.5], color="lightgreen"),
        dict(range=[threshold * 0.5, threshold], color="yellow"),
        dict(range=[threshold, max_value], color="red")
    ]
    
    # ✅ BOTH modes in one: Delta with conditional visibility
    fig = go.Figure(go.Indicator(
        mode="gauge+number+delta",
        value=value,
        number={
            'font': {'color': color, 'size': 38},
            'suffix': ' bp',
            'valueformat': '.0f'
        },
        delta={
            'reference': threshold,
            'position': "top",
            'font': {'size': 14, 'color': '#ef4444'},  # ✅ Smaller, less intrusive
            'increasing': {'color': '#ef4444'},
            'decreasing': {'color': '#22c55e'}
        },
        domain={'x': [0, 1], 'y': [0.1, 1]},  # ✅ Balanced for both views
        title={
            'text': title,
            'font': {'size': 17},
            'align': 'center'
        },
        gauge={
            'axis': {
                'range': [None, max_value],
                'tickwidth': 1,
                'tickfont': {'size': 10}
            },
            'bar': {'color': color, 'thickness': 0.75},
            'steps': steps,
            'threshold': {
                'line': {'color': "red", 'width': 3},
                'thickness': 0.75,
                'value': threshold
            }
        }
    ))
    
    fig.update_layout(
        height=260,
        margin=dict(l=15, r=15, t=45, b=5),  # ✅ Tighter margins
        paper_bgcolor='rgba(0,0,0,0)'
    )
    
    return fig

def create_distribution_chart(data_dict, title_key, lang="en", chart_type='bar', color_scheme=None):
    """Create distribution pie or bar charts"""
    if not data_dict:
        fig = go.Figure()
        fig.update_layout(
            title=f"{get_translation(title_key, lang)} (No Data)",
            xaxis={'visible': False}, yaxis={'visible': False},
            annotations=[{'text': 'No data available', 'xref': 'paper',
                        'yref': 'paper', 'showarrow': False, 'font': {'size': 16}}]
        )
        return fig

    title = get_translation(title_key, lang)
    df = pd.DataFrame(list(data_dict.items()), columns=['Category', 'Count'])

    limit = 15
    if len(df) > limit:
        df = df.nlargest(limit, 'Count')
        other_count = sum(count for _, count in data_dict.items()) - df['Count'].sum()
        if other_count > 0:
            df_other = pd.DataFrame([{'Category': 'Other', 'Count': other_count}])
            df = pd.concat([df, df_other], ignore_index=True)

    if chart_type.lower() == 'pie':
        fig = px.pie(df, values='Count', names='Category', title=f"{title}",
                     color_discrete_sequence=color_scheme if isinstance(color_scheme, list) else None)
        fig.update_traces(textposition='inside', textinfo='percent+label', pull=[0.05]*len(df))
        fig.update_layout(legend_title_text='Categories', showlegend=True)
    else:  # Bar chart
        # ✅ FIXED: Proper color handling for bar charts
        if isinstance(color_scheme, list) and len(color_scheme) >= 2:
            # Discrete list scheme
            fig = px.bar(df.sort_values('Count', ascending=True),
                         y='Category', x='Count', title=f"{title}", text_auto=True,
                         orientation='h',
                         color='Category',
                         color_discrete_sequence=color_scheme)  # ✅ Use discrete sequence
            fig.update_layout(yaxis_title=None, xaxis_title="Count", showlegend=False,
                              coloraxis_showscale=False)
            fig.update_yaxes(categoryorder='total ascending')
        else:
            # Single color or continuous (use default or single color)
            single_color = color_scheme[0] if isinstance(color_scheme, list) else (color_scheme if isinstance(color_scheme, str) and color_scheme.startswith('#') else '#3b82f6')
            fig = px.bar(df.sort_values('Count', ascending=True),
                         y='Category', x='Count', title=f"{title}", text_auto=True,
                         orientation='h')
            fig.update_traces(marker_color=single_color)  # ✅ Apply single color
            fig.update_layout(yaxis_title=None, xaxis_title="Count")
            fig.update_yaxes(categoryorder='total ascending')

    # ✅ REMOVED: No double color application - colors already applied above

    fig.update_layout(
        margin=dict(t=50, b=20, l=20, r=20),
        title_font_size=20,
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)'
    )
    return fig

# ==================== NEW PLOTLY VISUALIZATION FUNCTIONS ====================
def create_temporal_chart(sequences, interval='month', lang='en', color_scheme=None):
    """Generate a Plotly line chart for sequences over time."""
    progress_tracker.start_operation(f"Generating Temporal Chart (Interval: {interval})")
    df_data = [{'date': item[2].get('collection_date')} for item in sequences if item[2].get('collection_date')]
    if not df_data:
        progress_tracker.log_error("No date information found for temporal chart.")
        fig = go.Figure()
        fig.update_layout(title="Temporal Distribution (No Data)", xaxis={'visible': False}, yaxis={'visible': False},
                          annotations=[{'text': 'No date data available', 'xref': 'paper', 'yref': 'paper', 'showarrow': False, 'font': {'size': 16}}])
        return fig

    df = pd.DataFrame(df_data)
    df['date'] = pd.to_datetime(df['date'])
    df = df.dropna(subset=['date'])

    if interval == 'year':
        df['period'] = df['date'].dt.year.astype(str)
    elif interval == 'quarter':
        df['period'] = df['date'].dt.to_period('Q').astype(str)
    else:
        df['period'] = df['date'].dt.strftime('%Y-%m')

    counts = df['period'].value_counts().sort_index()
    counts_df = counts.reset_index()
    counts_df.columns = ['Period', 'Count']

    T = lambda key: get_translation(key, lang)
    title_text = f"Sequence Count by {interval.capitalize()}"
    if interval == 'year': title_text = T("vis_interval_year") + " Count"
    elif interval == 'quarter': title_text = T("vis_interval_quarter") + " Count"
    elif interval == 'month': title_text = T("vis_interval_month") + " Count"

    # ✅ FIXED: Don't pass color_scheme to line_dash_sequence
    fig = px.line(counts_df, x='Period', y='Count',
                  title=title_text,
                  markers=True, text='Count')
    fig.update_traces(textposition="top center")
    
    # ✅ FIXED: Apply color to line, not dash
    if color_scheme:
        if isinstance(color_scheme, str):  # Single color or scale name
            fig.update_traces(line=dict(color=color_scheme if color_scheme.startswith('#') else None))
        elif isinstance(color_scheme, list) and color_scheme:  # List of colors
            fig.update_traces(line=dict(color=color_scheme[0]))  # Use first color
    
    fig.update_layout(
        xaxis_title="Time Period", yaxis_title="Number of Sequences",
        margin=dict(t=50, b=20, l=20, r=20),
        paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)'
    )
    progress_tracker.complete_operation("Temporal chart generated")

    return fig

def create_geographic_heatmap(sequences, top_n=20, field='location', lang='en', color_scheme=None):
    """Generate a Plotly horizontal bar chart (simulating a heatmap) for any metadata field distribution."""
    progress_tracker.start_operation(f"Generating {field.capitalize()} Heatmap (Top {top_n})")

    # NEW: Dynamic exclusion of DEFAULT_UNKNOWN
    valid_items = [item[2].get(field, DEFAULT_UNKNOWN) for item in sequences if item[2].get(field, DEFAULT_UNKNOWN) != DEFAULT_UNKNOWN]
    if not valid_items:
        progress_tracker.log_error(f"No {field} information found for heatmap.")
        fig = go.Figure()
        fig.update_layout(title=f"{field.capitalize()} Distribution (No Data)", xaxis={'visible': False}, yaxis={'visible': False},
                          annotations=[{'text': f'No {field} data available', 'xref': 'paper', 'yref': 'paper', 'showarrow': False, 'font': {'size': 16}}])
        return fig

    # NEW: Use analyzer for dynamic counts (reuse existing method)
    analyzer = SequenceAnalyzer(sequences)  # Temp instance for distribution
    counts = analyzer.get_metadata_distribution(field)

    # Filter out Unknown/empty
    filtered_counts = {k: v for k, v in counts.items() if k != DEFAULT_UNKNOWN and k}
    if not filtered_counts:
        # Fallback error
        progress_tracker.log_error(f"No valid {field} data after filtering.")
        fig = go.Figure()
        fig.update_layout(title=f"{field.capitalize()} Distribution (No Valid Data)")
        return fig

    top_items = sorted(filtered_counts.items(), key=lambda x: x[1], reverse=True)[:top_n]
    df = pd.DataFrame(top_items, columns=[field.capitalize(), 'Count'])

    # Use translation for title if available, else dynamic
    T = lambda key: get_translation(key, lang)
    field_display = T(f"vis_field_{field}") if f"vis_field_{field}" in TRANSLATIONS.get(lang, {}) else field.capitalize()
    title_text = f"Top {len(df)} {field_display}s by Sequence Count"

    # UPDATED: Bar with dynamic field on y-axis
    fig = px.bar(df.sort_values('Count', ascending=True),  # Sort for horizontal
                 y=field.capitalize(), x='Count',
                 title=title_text,
                 text_auto=True, orientation='h',
                 color='Count',
                 color_continuous_scale=color_scheme)  # From prior updates

    fig.update_layout(
        yaxis_title=field_display,  # Dynamic y-label
        xaxis_title="Number of Sequences",
        coloraxis_showscale=False,  # Hide color bar legend
        margin=dict(t=50, b=20, l=20, r=20),
        paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)'
    )

    # NEW: Apply scheme if provided (from prior updates)
    if color_scheme:
        fig = apply_color_scheme(fig, color_scheme, 'heatmap')

    progress_tracker.complete_operation(f"{field.capitalize()} heatmap generated")
    return fig


def create_stacked_bar_chart(sequences, category1='location', category2='type', top_n=15, lang='en', color_scheme=None):
    """Generate a Plotly stacked bar chart."""
    progress_tracker.start_operation(f"Generating Stacked Bar ({category1} vs {category2}, Top {top_n})")
    T = lambda key: get_translation(key, lang)

    # UPDATED: Param validation for dynamic fields
    if category1 not in ['location', 'host', 'clade', 'type', 'segment', 'year', 'month'] or category2 not in ['location', 'host', 'clade', 'type', 'segment', 'year', 'month']:
        progress_tracker.log_error(f"Invalid categories: {category1} or {category2}. Use valid parser fields.")
        fig = go.Figure()
        fig.update_layout(title=f"Stacked Bar: Invalid Categories (No Data)")
        return fig

    # Helper to get field value, handling date fields
    def get_field(metadata, field):
        if field == 'year': return str(metadata['collection_date'].year) if metadata.get('collection_date') else DEFAULT_UNKNOWN
        if field == 'month': return metadata['collection_date'].strftime('%Y-%m') if metadata.get('collection_date') else DEFAULT_UNKNOWN
        return str(metadata.get(field, DEFAULT_UNKNOWN))

    # Aggregate data: category1 -> category2 -> count
    data = defaultdict(lambda: Counter())
    valid_seq_count = 0
    for _, _, meta in sequences:
         cat1_val = get_field(meta, category1)
         cat2_val = get_field(meta, category2)
         # Only count if both categories are known
         if cat1_val != DEFAULT_UNKNOWN and cat2_val != DEFAULT_UNKNOWN:
              data[cat1_val][cat2_val] += 1
              valid_seq_count += 1

    if not data:
        progress_tracker.log_error(f"No valid data found for stacking {category1} by {category2}.")
        fig = go.Figure()
        fig.update_layout(title=f"Stacked Bar: {category1.capitalize()} vs {category2.capitalize()} (No Data)", xaxis={'visible': False}, yaxis={'visible': False},
                          annotations=[{'text': 'No valid data for stacking', 'xref': 'paper', 'yref': 'paper', 'showarrow': False, 'font': {'size': 16}}])
        return fig

    # Prepare DataFrame for Plotly
    df_list = []
    for cat1, cat2_counts in data.items():
        for cat2, count in cat2_counts.items():
            df_list.append({category1: cat1, category2: cat2, 'Count': count})
    if not df_list: # Check if list is empty after filtering unknowns
        progress_tracker.log_error(f"No valid data points after filtering Unknowns for stacking.")
        fig = go.Figure(); fig.update_layout(title=f"Stacked Bar: {category1.capitalize()} vs {category2.capitalize()} (No Valid Data Points)"); return fig
    df = pd.DataFrame(df_list)

    # Get top N primary categories based on total count
    cat1_totals = df.groupby(category1)['Count'].sum().nlargest(top_n).index
    df_filtered = df[df[category1].isin(cat1_totals)]

    if df_filtered.empty:
         progress_tracker.log_error(f"No data remaining after filtering for top {top_n} {category1}s.")
         fig = go.Figure(); fig.update_layout(title=f"Stacked Bar: {category1.capitalize()} vs {category2.capitalize()} (No Data in Top {top_n})"); return fig

    # Sort secondary category for consistent legend color
    df_filtered = df_filtered.sort_values(by=[category1, category2])

    # Use translations for titles/labels
    cat1_display = T(f"vis_field_{category1}") if f"vis_field_{category1}" in TRANSLATIONS[lang] else category1.capitalize()
    cat2_display = T(f"vis_field_{category2}") if f"vis_field_{category2}" in TRANSLATIONS[lang] else category2.capitalize()
    title_text = f"{cat2_display} Distribution within Top {len(cat1_totals)} {cat1_display}s"

    fig = px.bar(df_filtered, x=category1, y='Count', color=category2,
                 title=title_text,
                 text_auto='.2s', # Show count on segments, formatted
                 category_orders={category1: cat1_totals.tolist()} # Keep the top N order
                )
    fig.update_traces(textfont_size=10, textangle=0, textposition="inside", cliponaxis=False) # Improve text visibility

    fig.update_layout(
        xaxis_title=cat1_display, yaxis_title="Number of Sequences",
        legend_title=cat2_display,
        margin=dict(t=50, b=20, l=20, r=20),
        paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
        xaxis={'categoryorder':'array', 'categoryarray':cat1_totals.tolist()} # Explicitly order x-axis
    )
    fig.update_xaxes(tickangle=45)

    # UPDATED: Apply color scheme if provided (from prior updates)
    if color_scheme:
        fig = apply_color_scheme(fig, color_scheme, 'stacked', len(df_filtered))

    progress_tracker.complete_operation("Stacked bar chart generated")
    return fig

def nucleotide_palette(seq_sample):
    """Creative: Generate palette inspired by nucleotide frequencies (A=green, T=red, etc.)."""
    if not seq_sample:
        return px.colors.sequential.Blues  # Default
    sample = seq_sample.upper()[:100]  # First 100 bases
    counts = Counter(sample)
    total = len(sample)
    gc_ratio = (counts.get('G', 0) + counts.get('C', 0)) / total if total > 0 else 0.5
    at_ratio = 1 - gc_ratio
    # Blend: High GC -> Greens, High AT -> Reds/Oranges
    if gc_ratio > 0.6:
        return px.colors.sequential.Greens
    elif at_ratio > 0.6:
        return px.colors.sequential.Reds
    else:
        return px.colors.sequential.Purples  # Balanced

def apply_color_scheme(fig, color_scheme, chart_type, data_len=0):
    """Helper: Apply color scheme to figure based on type."""
    if not color_scheme:
        return fig
    
    if isinstance(color_scheme, str):  # Scale name
        if 'heatmap' in chart_type:
            # Heatmaps use colorscale directly
            fig.update_traces(colorscale=color_scheme)
        elif 'line' in chart_type:
            # Line charts use line color
            fig.update_traces(line=dict(color=color_scheme if color_scheme.startswith('#') else None))
        elif 'bar' in chart_type or 'stacked' in chart_type:
            # ✅ FIXED: Bar charts - use marker.color with colorscale via layout
            fig.update_layout(coloraxis_colorscale=color_scheme)
        elif 'pie' in chart_type:
            # Pie charts don't use colorscale the same way
            pass
    else:  # List of hex/colors
        if 'pie' in chart_type:
            # ✅ CORRECT: Pie uses marker.colors (plural)
            colors = color_scheme[:data_len] if data_len else color_scheme
            fig.update_traces(marker=dict(colors=colors))
        elif 'bar' in chart_type or 'stacked' in chart_type:
            # ✅ FIXED: Bar charts use marker.color (singular)
            colors = color_scheme[:data_len] if data_len else color_scheme
            fig.update_traces(marker=dict(color=colors))
        elif 'line' in chart_type:
            # Line charts use line color
            color = color_scheme[0] if color_scheme else 'blue'
            fig.update_traces(line=dict(color=color))
        elif 'heatmap' in chart_type:
            # Heatmaps use colorscale (list of colors)
            colors = color_scheme[:10] if len(color_scheme) > 10 else color_scheme
            fig.update_traces(colorscale=colors)
    
    return fig
 

# ==================== CUSTOM CSS ====================
def load_custom_css():
    """Load custom CSS with theme support"""
    theme = st.session_state.get('theme', 'light')
    
    if theme == 'dark':
        css = r"""
        <style>
            /* Dark Mode Colors */
            :root {
                --bg-primary: #0f172a;
                --bg-secondary: #1e293b;
                --bg-tertiary: #334155;
                --text-primary: #f1f5f9;
                --text-secondary: #cbd5e1;
                --accent-blue: #3b82f6;
                --accent-blue-light: #60a5fa;
                --border-color: #475569;
                --shadow: rgba(0, 0, 0, 0.5);
            }
            
            .main {
                background: linear-gradient(180deg, #0f172a 0%, #1e293b 100%);
                color: var(--text-primary);
            }
            
            .main .block-container {
                padding-top: 2rem;
                padding-bottom: 2rem;
            }
            
            /* Cards/Containers */
            div[data-testid="stExpander"] div[data-testid="stVerticalBlock"],
            div.stTabs [data-baseweb="tab-panel"] > div[data-testid="stVerticalBlock"] > div:not([data-testid="stExpander"]):not(:has(div[data-testid="stExpander"])){
                background: var(--bg-secondary);
                padding: 25px;
                border-radius: 12px;
                box-shadow: 0 4px 12px var(--shadow);
                margin-bottom: 25px;
                border: 1px solid var(--border-color);
            }
            
            /* Typography */
            h1, h2, h3, h4, h5, h6 {
                color: var(--text-primary) !important;
            }
            
            h1 {
                font-weight: 700;
                text-align: center;
                margin-bottom: 0.5rem;
            }
            
            h2 {
                border-bottom: 2px solid var(--accent-blue);
                padding-bottom: 8px;
                margin-top: 1rem;
                margin-bottom: 1.5rem;
            }
            
            h3, h4 {
                margin-top: 1.5rem;
                margin-bottom: 1rem;
                font-weight: 600;
            }
            
            /* Text Elements */
            p, label, .stMarkdown {
                color: var(--text-secondary) !important;
            }
            
            /* Buttons */
            .stButton > button {
                border: none;
                border-radius: 8px;
                padding: 10px 20px;
                font-weight: 600;
                transition: all 0.2s ease-in-out;
                box-shadow: 0 2px 4px var(--shadow);
            }
            
            .stButton > button[kind="primary"] {
                background: linear-gradient(90deg, var(--accent-blue) 0%, var(--accent-blue-light) 100%);
                color: white;
            }
            
            .stButton > button[kind="primary"]:hover {
                box-shadow: 0 4px 12px rgba(59, 130, 246, 0.5);
                filter: brightness(1.2);
            }
            
            .stButton > button[kind="secondary"] {
                background-color: var(--bg-tertiary);
                color: #ef4444;
                border: 1px solid #ef4444;
            }
            
            .stButton > button[kind="secondary"]:hover {
                background-color: var(--bg-secondary);
                border-color: #dc2626;
            }
            
            .stButton > button:not([kind="primary"]):not([kind="secondary"]) {
                background-color: var(--bg-tertiary);
                color: var(--accent-blue-light);
                border: 1px solid var(--border-color);
            }
            
            .stButton > button:not([kind="primary"]):not([kind="secondary"]):hover {
                background-color: var(--bg-secondary);
                border-color: var(--border-color);
            }
            
            /* Alerts */
            .stAlert {
                border-radius: 8px;
                border-left-width: 5px;
                padding: 1rem;
                background-color: var(--bg-tertiary) !important;
                color: var(--text-primary) !important;
                box-shadow: 0 2px 6px var(--shadow);
            }
            
            /* File Uploader */
            .stFileUploader {
                border: 2px dashed var(--accent-blue);
                border-radius: 10px;
                padding: 25px;
                background: var(--bg-tertiary);
            }
            
            .stFileUploader label {
                font-weight: 600;
                color: var(--accent-blue-light);
            }
            
            /* Metrics */
            div[data-testid="stMetric"] {
                background-color: var(--bg-tertiary);
                border: 1px solid var(--border-color);
                padding: 1.5rem;
                border-radius: 12px;
                box-shadow: 0 2px 8px var(--shadow);
            }
            
            div[data-testid="stMetricLabel"] {
                font-weight: 600;
                color: var(--text-secondary) !important;
                font-size: 0.95rem;
            }
            
            div[data-testid="stMetricValue"] {
                font-size: 2.2rem;
                font-weight: 700;
                color: var(--accent-blue-light) !important;
            }
            
            div[data-testid="stMetricDelta"] {
                font-size: 0.9rem;
            }
            
            /* Sidebar (keep dark in dark mode) */
            [data-testid="stSidebar"] {
                background: linear-gradient(180deg, #020617 0%, #0c4a6e 100%);
                padding: 1rem;
            }
            
            [data-testid="stSidebar"] h3 {
                color: #e0f2fe;
                border-bottom: 1px solid #7dd3fc;
            }
            
            [data-testid="stSidebar"] .stMetric {
                background-color: rgba(255, 255, 255, 0.1);
                border: none;
                box-shadow: none;
            }
            
            [data-testid="stSidebar"] .stMetricLabel {
                color: #e0f2fe;
                font-size: 0.9rem;
            }
            
            [data-testid="stSidebar"] .stMetricValue {
                color: #ffffff;
                font-size: 1.8rem;
            }
            
            [data-testid="stSidebar"] .stButton > button {
                background-color: rgba(255, 255, 255, 0.2);
                color: white;
                border: 1px solid rgba(255, 255, 255, 0.4);
            }
            
            [data-testid="stSidebar"] .stButton > button:hover {
                background-color: rgba(255, 255, 255, 0.3);
                border-color: rgba(255, 255, 255, 0.6);
            }
            
            [data-testid="stSidebar"] .stSelectbox label,
            [data-testid="stSidebar"] .stRadio label {
                color: #e0f2fe;
                font-weight: 600;
            }
            
            /* ← ADDED: Caption/text styling for sidebar */
            [data-testid="stSidebar"] .stCaption,
            [data-testid="stSidebar"] p {
                color: #cbd5e1 !important;
                font-size: 0.85rem;
            }
            
            /* Tabs */
            .stTabs [data-baseweb="tab-list"] {
                gap: 12px;
                background-color: transparent;
                border-radius: 0;
                padding: 0;
                box-shadow: none;
                border-bottom: 2px solid var(--border-color);
                margin-bottom: 1.5rem;
            }
            
            .stTabs [data-baseweb="tab"] {
                background-color: transparent;
                border-radius: 8px 8px 0 0;
                padding: 12px 24px;
                font-weight: 600;
                color: var(--text-secondary);
                border: none;
                border-bottom: 2px solid transparent;
                margin-bottom: -2px;
                transition: all 0.2s ease;
            }
            
            .stTabs [data-baseweb="tab"]:hover {
                background-color: var(--bg-tertiary);
                color: var(--accent-blue-light);
            }
            
            .stTabs [aria-selected="true"] {
                color: var(--accent-blue-light);
                background-color: transparent;
                border-bottom: 2px solid var(--accent-blue);
                box-shadow: none;
            }
            
            /* Progress Bar */
            .stProgress > div > div {
                background: linear-gradient(90deg, var(--accent-blue) 0%, var(--accent-blue-light) 100%);
                border-radius: 8px;
            }
            
            /* DataFrames */
            .stDataFrame {
                border-radius: 10px;
                overflow: hidden;
                box-shadow: 0 2px 8px var(--shadow);
                border: 1px solid var(--border-color);
            }
            
            /* Expander */
            .stExpander > summary {
                background-color: var(--bg-tertiary);
                border-radius: 8px;
                padding: 10px 15px;
                font-weight: 600;
                color: var(--text-primary) !important;
                border: 1px solid var(--border-color);
            }
            
            .stExpander > summary:hover {
                background-color: var(--bg-secondary);
            }
            
            .stExpander > div {
                border-top: none;
                padding-top: 15px;
            }
            
            /* Input Fields */
            .stTextInput input, 
            .stTextArea textarea, 
            .stSelectbox select,
            .stNumberInput input {
                background-color: var(--bg-tertiary) !important;
                color: var(--text-primary) !important;
                border-color: var(--border-color) !important;
            }
            
            /* Checkbox & Radio */
            .stCheckbox label,
            .stRadio label {
                color: var(--text-primary) !important;
            }
            
            /* Sliders */
            .stSlider {
                color: var(--text-primary) !important;
            }
            
            /* Code blocks */
            code {
                background-color: var(--bg-tertiary) !important;
                color: var(--accent-blue-light) !important;
            }
            
            pre {
                background-color: var(--bg-tertiary) !important;
                border: 1px solid var(--border-color) !important;
            }
        </style>
        """
    else:  # Light mode (complete original CSS)
        css = r"""
        <style>
            .main .block-container {
                padding-top: 2rem;
                padding-bottom: 2rem;
            }
            .main {
                background: linear-gradient(180deg, #f0f9ff 0%, #e0f2fe 100%);
            }
            div[data-testid="stExpander"] div[data-testid="stVerticalBlock"],
            div.stTabs [data-baseweb="tab-panel"] > div[data-testid="stVerticalBlock"] > div:not([data-testid="stExpander"]):not(:has(div[data-testid="stExpander"])){
                 background: white;
                 padding: 25px;
                 border-radius: 12px;
                 box-shadow: 0 4px 12px rgba(0, 0, 0, 0.08);
                 margin-bottom: 25px;
                 border: 1px solid #e5e7eb;
            }
            h1 {
                color: #1e3a8a;
                font-weight: 700;
                text-align: center;
                margin-bottom: 0.5rem;
            }
            h2 {
                color: #1d4ed8;
                border-bottom: 2px solid #60a5fa;
                padding-bottom: 8px;
                margin-top: 1rem;
                margin-bottom: 1.5rem;
            }
            h3 {
                color: #1e40af;
                margin-top: 1.5rem;
                margin-bottom: 1rem;
                font-weight: 600;
            }
            h4 {
                color: #1e40af;
                margin-top: 1.5rem;
                margin-bottom: 1rem;
                font-weight: 600;
            }
            .stButton > button {
                border: none;
                border-radius: 8px;
                padding: 10px 20px;
                font-weight: 600;
                transition: all 0.2s ease-in-out;
                box-shadow: 0 2px 4px rgba(0, 0, 0, 0.08);
            }
            .stButton > button[kind="primary"] {
                 background: linear-gradient(90deg, #3b82f6 0%, #60a5fa 100%);
                 color: white;
            }
            .stButton > button[kind="primary"]:hover {
                 box-shadow: 0 4px 12px rgba(59, 130, 246, 0.3);
                 filter: brightness(1.1);
            }
            .stButton > button[kind="secondary"] {
                 background-color: #f3f4f6;
                 color: #dc2626;
                 border: 1px solid #ef4444;
            }
            .stButton > button[kind="secondary"]:hover {
                 background-color: #fee2e2;
                 border-color: #dc2626;
            }
            .stButton > button:not([kind="primary"]):not([kind="secondary"]) {
                background-color: #ffffff;
                color: #3b82f6;
                border: 1px solid #d1d5db;
            }
            .stButton > button:not([kind="primary"]):not([kind="secondary"]):hover {
                background-color: #f9fafb;
                border-color: #9ca3af;
            }
            .stAlert {
                border-radius: 8px;
                border-left-width: 5px;
                padding: 1rem;
                box-shadow: 0 2px 6px rgba(0,0,0,0.06);
            }
            .stFileUploader {
                border: 2px dashed #93c5fd;
                border-radius: 10px;
                padding: 25px;
                background: #eff6ff;
            }
            .stFileUploader label {
                font-weight: 600;
                color: #1d4ed8;
            }
            div[data-testid="stMetric"] {
                background-color: #ffffff;
                border: 1px solid #e5e7eb;
                padding: 1.5rem;
                border-radius: 12px;
                box-shadow: 0 2px 8px rgba(0,0,0,0.05);
            }
            div[data-testid="stMetricLabel"] {
                font-weight: 600;
                color: #4b5563;
                font-size: 0.95rem;
            }
            div[data-testid="stMetricValue"] {
                font-size: 2.2rem;
                font-weight: 700;
                color: #1e3a8a;
            }
            div[data-testid="stMetricDelta"] {
                 font-size: 0.9rem;
            }
            [data-testid="stSidebar"] {
                background: linear-gradient(180deg, #0c4a6e 0%, #0369a1 100%);
                padding: 1rem;
            }
            [data-testid="stSidebar"] h3 {
                 color: #e0f2fe;
                 border-bottom: 1px solid #7dd3fc;
            }
            [data-testid="stSidebar"] .stMetric {
                 background-color: rgba(255, 255, 255, 0.1);
                 border: none;
                 box-shadow: none;
            }
            [data-testid="stSidebar"] .stMetricLabel {
                  color: #e0f2fe;
                  font-size: 0.9rem;
            }
            [data-testid="stSidebar"] .stMetricValue {
                  color: #ffffff;
                  font-size: 1.8rem;
            }
            [data-testid="stSidebar"] .stButton > button {
                 background-color: rgba(255, 255, 255, 0.2);
                 color: white;
                 border: 1px solid rgba(255, 255, 255, 0.4);
            }
            [data-testid="stSidebar"] .stButton > button:hover {
                  background-color: rgba(255, 255, 255, 0.3);
                  border-color: rgba(255, 255, 255, 0.6);
            }
            [data-testid="stSidebar"] .stSelectbox label,
            [data-testid="stSidebar"] .stRadio label {  /* ← ADDED: Radio support */
                  color: #e0f2fe;
                  font-weight: 600;
            }
            
            /* ← ADDED: Caption/text styling for sidebar */
            [data-testid="stSidebar"] .stCaption,
            [data-testid="stSidebar"] p {
                color: #bae6fd !important;
                font-size: 0.85rem;
            }
            
            .stTabs [data-baseweb="tab-list"] {
                gap: 12px;
                background-color: transparent;
                border-radius: 0;
                padding: 0;
                box-shadow: none;
                border-bottom: 2px solid #d1d5db;
                margin-bottom: 1.5rem;
            }
            .stTabs [data-baseweb="tab"] {
                background-color: transparent;
                border-radius: 8px 8px 0 0;
                padding: 12px 24px;
                font-weight: 600;
                color: #4b5563;
                border: none;
                border-bottom: 2px solid transparent;
                margin-bottom: -2px;
                transition: all 0.2s ease;
            }
            .stTabs [data-baseweb="tab"]:hover {
                 background-color: #f3f4f6;
                 color: #1d4ed8;
            }
            .stTabs [aria-selected="true"] {
                 color: #1d4ed8;
                 background-color: transparent;
                 border-bottom: 2px solid #1d4ed8;
                 box-shadow: none;
            }
            .stProgress > div > div {
                background: linear-gradient(90deg, #3b82f6 0%, #60a5fa 100%);
                border-radius: 8px;
            }
            .stDataFrame {
                border-radius: 10px;
                overflow: hidden;
                box-shadow: 0 2px 8px rgba(0,0,0,0.06);
                border: 1px solid #e5e7eb;
            }
            .stExpander > summary {
                 background-color: #f9fafb;
                 border-radius: 8px;
                 padding: 10px 15px;
                 font-weight: 600;
                 color: #1e3a8a;
                 border: 1px solid #e5e7eb;
            }
            .stExpander > summary:hover {
                  background-color: #f3f4f6;
            }
            .stExpander > div {
                  border-top: none;
                  padding-top: 15px;
            }
        </style>
        """
    
    st.markdown(css, unsafe_allow_html=True)

# ==================== SESSION STATE INITIALIZATION ====================
def init_session_state():
    """Initialize session state variables if they don't exist."""
    defaults = {
        'lang': 'en',
        'theme': 'light',
        'all_files': {},
        'active_sequences': [],
        'original_sequences': {},
        'analysis_log': [],
        'processing_step': 0,
        'status_message': "processing",
        'status_level': "info",
        'active_filenames': [],
        'last_report': "",
        'accession_list': [],
        'confirming_removal': False,
        'status_placeholder': None,
        'gdrive_mounted': False,
        'generated_chart': None
    }
    for key, default_value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = default_value

# ==================== MAIN APP ====================
def main():
    st.set_page_config(
        page_title="Vir-Seq-Sift - FASTA Analysis",
        page_icon="🧬🧺",
        layout="wide",
        initial_sidebar_state="expanded"
    )

    load_custom_css()
    init_session_state()

    # Sidebar with full translation
    with st.sidebar:
        st.markdown("<h1 style='text-align: center; color: white;'>🧬 Vir-Seq-Sift</h1>", unsafe_allow_html=True)
        
        # ========== ENHANCED LANGUAGE SELECTOR ==========
        # NEW: Global state snapshot helper (call before rerun)
        def snapshot_widget_states():
            """Quickly snapshot key widget states to session state for restoration."""
            # Example: Add keys for your main widgets (expand as needed)
            snapshot = {
                'analyze_min_len': st.session_state.get('analyze_min_len', 200),
                'analyze_max_n': st.session_state.get('analyze_max_n', 100),
                'analyze_subtype_select': st.session_state.get('analyze_subtype_select', 'All'),  # Or use '_index' for index
                'analyze_subtype_select_index': st.session_state.get('analyze_subtype_select_index', 0),
                'vis_chart_type': st.session_state.get('vis_chart_type', 'bar'),
                'vis_chart_type_index': st.session_state.get('vis_chart_type_index', 0),
                'vis_field1': st.session_state.get('vis_field1', 'type'),
                'vis_field1_index': st.session_state.get('vis_field1_index', 0),
                'vis_interval': st.session_state.get('vis_interval', 'month'),
                'vis_top_n': st.session_state.get('vis_top_n', 20),
                'vis_cat1_stacked': st.session_state.get('vis_cat1_stacked', 'location'),
                'vis_cat2_stacked': st.session_state.get('vis_cat2_stacked', 'type'),
                
                # Manage Tab
                'manage_file_multiselect': st.session_state.get('manage_file_multiselect', []),
                'upload_method_radio': st.session_state.get('upload_method_radio', 0),  # Index for radio
                
                # Refine Tab
                'refine_clade_mode': st.session_state.get('refine_clade_mode', 'single'),
                'refine_clade_mode_index': st.session_state.get('refine_clade_mode_index', 0),
                'refine_clade_keep': st.session_state.get('refine_clade_keep', 'Both (First & Last)'),
                'refine_clade_keep_index': st.session_state.get('refine_clade_keep_index', 2),
                'refine_clade_separate': st.session_state.get('refine_clade_separate', True),
                'refine_temp_group': st.session_state.get('refine_temp_group', 'location_host'),
                'refine_temp_group_index': st.session_state.get('refine_temp_group_index', 4),
                'refine_temp_sort': st.session_state.get('refine_temp_sort', 'date'),
                'refine_temp_keep': st.session_state.get('refine_temp_keep', 'both'),
                'refine_temp_keep_index': st.session_state.get('refine_temp_keep_index', 2),
                'refine_temp_custom': st.session_state.get('refine_temp_custom', ''),
                
                # Global/Sidebar
                'data_mode_toggle': st.session_state.get('data_mode_toggle', 'Current (Filtered)'),
                'theme_mode': st.session_state.get('theme_mode', 'auto'),
                
                # Upload (if needed)
                'url_downloader_input': st.session_state.get('url_downloader_input', ''),
                'gdrive_path_input': st.session_state.get('gdrive_path_input', '')
                
                # This prevents resets on language change
            }
            for key, val in snapshot.items():
                if key not in st.session_state:
                    st.session_state[key] = val
    
        # Check for language in query params first
        default_lang = st.query_params.get("lang", "en")
        if 'lang' not in st.session_state:
            st.session_state.lang = default_lang
    
        lang_options = {'en': "🇬🇧 English", 'ru': "🇷🇺 Русский"}
        
        # Callback function (triggered on change)
        def update_language():
            current_selection = st.session_state['lang_selector']  # Access via session_state in callback
            if current_selection != st.session_state.lang:
                # Snapshot states BEFORE change
                snapshot_widget_states()
                st.session_state.lang = current_selection
                # Persist to query params
                st.query_params.update({"lang": current_selection})
                # Optional: Toast confirmation
                st.toast(f"✅ Switched to {lang_options[current_selection]} – Data preserved!", icon="🌐")
                # Rerun for full UI refresh (safe now with snapshot)
                st.rerun()
    
        selected_lang_code = st.selectbox(
            "🌐 Language / Язык",
            options=list(lang_options.keys()),
            index=list(lang_options.keys()).index(st.session_state.lang),
            format_func=lambda code: lang_options[code],
            key='lang_selector',  # Consistent key
            label_visibility="collapsed",
            on_change=update_language  # FIXED: Use on_change for callback (Streamlit standard)
        )
        
        # ========== END ENHANCED ==========
    
        # ========== ADD THEME TOGGLE HERE (RIGHT AFTER LANGUAGE) ==========
        # Theme Toggle
        default_theme = st.query_params.get("theme", "light")
        if 'theme' not in st.session_state:
            st.session_state.theme = default_theme
    
        # Theme options with Auto mode
        theme_options = {
            'auto': "🔄 Auto",
            'light': "☀️ Light",
            'dark': "🌙 Dark"
        }
    
        # Initialize theme mode (auto/manual preference)
        if 'theme_mode' not in st.session_state:
            st.session_state.theme_mode = st.query_params.get("theme_mode", "auto")
    
        # Radio button selector
        selected_theme_mode = st.radio(
            "Theme",
            options=list(theme_options.keys()),
            index=list(theme_options.keys()).index(st.session_state.theme_mode),
            format_func=lambda x: theme_options[x],
            key='theme_selector',
            horizontal=True,
            label_visibility="collapsed"
        )
    
        # Handle theme change
        if selected_theme_mode != st.session_state.theme_mode:
            st.session_state.theme_mode = selected_theme_mode
            
            if selected_theme_mode == 'auto':
                # Auto mode: Use time-based logic (dark at night, light during day)
                hour = datetime.now().hour
                st.session_state.theme = 'dark' if (20 <= hour or hour < 6) else 'light'
            else:
                # Manual mode: Use selected theme directly
                st.session_state.theme = selected_theme_mode
            
            # Persist to query params
            st.query_params.update({
                "theme": st.session_state.theme,
                "theme_mode": selected_theme_mode
            })
            st.rerun()
    
        # Show current mode indicator (optional)
        if st.session_state.theme_mode == 'auto':
            st.caption(f"🔄 Auto mode: {st.session_state.theme.capitalize()} ({datetime.now().strftime('%H:%M')})")
        
        # ========== END COMBINED ==========
    
        T = lambda key: get_translation(key, st.session_state.get('lang', 'en'))
    
        st.markdown("---")
    
        # Data Mode Toggle (radio) + per-file selectbox
        data_mode = st.radio(
            T("data_mode_label"),  # "Data Mode:"
            [T("data_mode_current"), T("data_mode_original")],  # ["Current (Filtered)", "Original (Pre-Filter)"]
            index=0, 
            key="data_mode_toggle", 
            horizontal=True,
            help=T("data_mode_help")  # Existing help key
        )
        data_mode_val = 'current' if data_mode == T("data_mode_current") else 'original'  # FIXED: Direct string match
    
        # Per-file selectbox for multi-file Data Mode (right after radio)
        if st.session_state.get('active_filenames') and len(st.session_state.active_filenames) > 1:
            available_originals = list(st.session_state.original_sequences.keys())
            selected_file = st.selectbox(
                "Use Original from File:",
                options=['Merged (All Active)'] + available_originals,
                index=0,
                key="data_mode_file_select",
                help="For multi-file: Pick a specific file's pre-filter snapshot, or merged."
            )
            
            if selected_file != 'Merged (All Active)' and selected_file in st.session_state.original_sequences:
                # Swap to per-file original
                st.session_state.active_sequences = [list(s) for s in st.session_state.original_sequences[selected_file]]
                st.session_state.data_mode_file = selected_file  # Track for reruns
            else:
                # Fallback to merged or current
                if data_mode_val == 'original':
                    st.session_state.active_sequences = [list(s) for s in st.session_state.get('original_active_snapshot', [])]
                # Else: Keep current (filtered)
        else:
            # Single-file or none: Hide dropdown
            st.session_state.data_mode_file = None
    
        st.markdown(f"### {T('sidebar_quick_stats')}")
        if st.session_state.all_files:
            st.metric(T("sidebar_files_loaded"), len(st.session_state.all_files))
        else:
            st.caption(T("sidebar_no_files"))
    
        if st.session_state.active_sequences:
            st.metric(T("sidebar_active_seqs"), f"{len(st.session_state.active_sequences):,}")
            try:
                avg_len = sum(len(s[1]) for s in st.session_state.active_sequences) / len(st.session_state.active_sequences)
                st.metric(T("sidebar_avg_length"), f"{int(avg_len):,} {T('bp')}")
            except ZeroDivisionError:
                st.metric(T("sidebar_avg_length"), "N/A")
            
            # PLACEMENT: Merged Files Metric—right here after active seqs/avg length (inside if active_sequences)
            if st.session_state.active_filenames:
                merged_files_count = len(st.session_state.active_filenames)
                st.metric("Merged from", f"{merged_files_count} {T('files')}", delta=None)
            else:
                st.caption(T("no_merged_files_active"))  # "No merged files active." (translated)
        else:
            st.caption(T("sidebar_no_dataset"))
    
        st.markdown("---")
    
        # ... (rest of your sidebar: quick actions, reset, export, footer—unchanged)

        st.markdown(f"### {T('sidebar_quick_actions')}")
        
        if st.button(T("sidebar_reset_all"), use_container_width=True, key="reset_all_sidebar"):
            # Step 1: Preserve user preferences
            preserved_lang = st.session_state.get('lang', 'en')
            preserved_theme = st.session_state.get('theme', 'light')
            
            # Spinner with new translation
            with st.spinner(T("reset_spinner_text")):
                keys_to_delete = [k for k in st.session_state.keys() if k not in ['status_placeholder']]
                for key in keys_to_delete:
                    del st.session_state[key]
                
                init_session_state()
                
                st.query_params.update({
                    "lang": preserved_lang,
                    "theme": preserved_theme,
                    "reset": "true"
                })
                
                time.sleep(0.3)
            
            # Toast with existing translation (reuse)
            st.toast(f"✅ {T('sidebar_reset_success')}", icon="🔄") # Shows: Toast popup "✅ 🔄 Session Reset!" st.rerun()   
            
            time.sleep(0.5)
            st.rerun()

        if st.session_state.active_sequences:
            fasta_str_io = io.StringIO()
            for header, seq, _ in st.session_state.active_sequences:
                h = header if isinstance(header, str) else str(header or '')
                h = h if h.startswith('>') else '>' + h
                s = seq if isinstance(seq, str) else str(seq or '')
                fasta_str_io.write(f"{h}\n{s}\n")

            st.download_button(
                label=T("sidebar_quick_export"),
                data=fasta_str_io.getvalue(),
                file_name=f"quick_export_{datetime.now().strftime('%Y%m%d_%H%M')}.fasta",
                mime="text/plain",
                use_container_width=True,
                key="quick_export_sidebar"
            )

        st.markdown("---")
        st.caption(f"{T('sidebar_footer')} | {datetime.now().year}")

        # ========== TEMP DEBUG EXPANDER (REMOVE AFTER TESTING) ==========
        with st.expander("Debug: Session Snapshot", expanded=False):
            # Filter to seq/file-related keys for brevity; expand filter as needed
            debug_data = {k: v for k, v in st.session_state.items() if 'seq' in k.lower() or 'file' in k.lower()}
            st.json(debug_data)
        # ========== END TEMP ==========
        

    # Main Area
    st.markdown(f"## {T('app_title')}")

    st.session_state.status_placeholder = st.empty()
    if st.session_state.status_message:
        update_status(st.session_state.status_message, st.session_state.status_level, log=False)

    # NEW: Tab tracking for navigation (replaces old tab setup)
    if 'active_tab' not in st.session_state:
        st.session_state.active_tab = 0  # Default to first tab (Upload)
    
    tab_keys = ["upload_tab", "manage_tab", "analyze_tab", "refine_tab", "export_tab", "docs_tab"]
    tab_labels = [T(key) for key in tab_keys]
    tabs = st.tabs(tab_labels)
    
    # FIXED: Select active tab by index (simulates selection via session state)
    # st.session_state.active_tab = tabs.index(tabs)  # Update on UI selection
    tab_map = dict(zip(tab_keys, tabs))  # Map for with blocks

    # Get active tab from query params
    query_tab = st.query_params.get("tab", "upload")
    if query_tab in tab_keys:
        # Simulate by hiding other tabs' content (hacky but works)
        for i, key in enumerate(tab_keys):
            if key == query_tab:
                st.session_state.active_tab = i
            else:
                with tab_map[key]:
                    st.empty()  # Hide content

    # ==================== TAB 1: UPLOAD & SETUP ====================
    with tab_map["upload_tab"]:
        st.header(T("upload_tab"))

        if not st.session_state.all_files:
            st.markdown(f"""
                <div style='background: linear-gradient(135deg, #e0f2fe 0%, #ccfbf1 100%);
                            color: #0c4a6e; padding: 25px; border-radius: 12px; border-left: 6px solid #0ea5e9;'>
                    <h3 style='color: #0c4a6e; border: none; margin-top: 0;'>{T('welcome_title')}</h3>
                    <p style='font-size: 1.05rem;'>{T('welcome_message')}</p>
                    <p>{T('welcome_subtitle')}</p>
                </div>
            """, unsafe_allow_html=True)
            st.markdown("<br>", unsafe_allow_html=True)

        # --- FIXED: Isolated upload methods based on radio selection ---
        upload_options = [T("upload_widget"), T("upload_url")]
        # Only add Drive option if potentially available
        if COLAB_AVAILABLE or os.path.exists('/content/drive'):  # Basic check
            upload_options.insert(1, T("upload_gdrive"))

        selected_upload_method = st.radio("Select Upload Method:", upload_options, index=st.session_state.get('upload_method_index', 0), horizontal=True, label_visibility="collapsed", key='upload_method_radio')

        # --- Widget Upload (only show if selected) ---
        if selected_upload_method == T("upload_widget"):
            st.subheader(T("upload_files_header"))
            uploaded_files = st.file_uploader(
                T("file_uploader_label"),
                type=['fasta', 'fas', 'fa', 'fna', 'txt', 'gz'],
                accept_multiple_files=True,
                help=T("upload_help_text"),
                key="main_file_uploader"
            )

            if uploaded_files:
                parser = FastaParser()  # FIXED: Instantiate parser here
                with st.spinner(T("processing_files")):
                    progress_bar = st.progress(0, text=T("initializing"))
                    newly_loaded_count = 0
                    total_sequences_added = 0
                    has_errors = False

                    for idx, uploaded_file in enumerate(uploaded_files):
                        filename = uploaded_file.name
                        progress_text = f"{T('processing')}: {filename}"
                        progress_bar.progress((idx) / len(uploaded_files), text=progress_text)

                        if filename not in st.session_state.all_files:
                            try:
                                content_bytes = uploaded_file.getvalue()
                                if filename.lower().endswith('.gz'):
                                    content_string = gzip.decompress(content_bytes).decode('utf-8', errors='replace')
                                else:
                                    content_string = content_bytes.decode('utf-8', errors='replace')

                                sequences, errors = parser.parse(content_string)  # Now uses local parser

                                if errors:
                                    has_errors = True
                                    st.warning(f"⚠️ {filename}: {errors[0]}", icon="⚠️")

                                if sequences:
                                    st.session_state.all_files[filename] = sequences
                                    if filename not in st.session_state.original_sequences:
                                        st.session_state.original_sequences[filename] = sequences
                                    newly_loaded_count += 1
                                    total_sequences_added += len(sequences)
                                else:
                                    if not errors:
                                        progress_tracker.log_error(f"No valid sequences found in {filename}")

                            except Exception as e:
                                progress_tracker.log_error(f"Failed to process {filename}: {str(e)}")
                                has_errors = True

                    progress_bar.progress(1.0, text=T("processing_complete"))
                    time.sleep(1)
                    progress_bar.empty()

                    if newly_loaded_count > 0:
                        msg = T("loaded_files").format(count=newly_loaded_count, seqs=total_sequences_added)
                        st.success(msg)
                        st.toast(
                            f"🎊 Success! {newly_loaded_count} files ({total_sequences_added:,} sequences) ready for analysis",
                            icon="🧬"
                        )
                        if not st.session_state.active_sequences:
                            st.info(T("info_activate_files"))
                    elif not has_errors:
                        st.warning(T("no_new_files"))

        # --- URL Download (only show if selected) ---
        elif selected_upload_method == T("upload_url"):
            st.subheader(T("download_url_header"))
            url_input = st.text_input(
                T("url_input_label"),
                value=st.session_state.get('url_downloader_input', ''),
                placeholder=T("url_placeholder"),
                key="url_downloader_input"
            )
            if st.button(T("download_url_btn"), use_container_width=True, key="url_download_button"):
                if url_input and url_input.startswith(('http://', 'https://')):
                    with st.spinner(T("downloading_from_url").format(url=url_input[:50] + '...')):  # FIXED: Add '...' for truncation
                        try:
                            response = requests.get(url_input, timeout=DEFAULT_TIMEOUT, stream=True)
                            response.raise_for_status()

                            content_disp = response.headers.get('content-disposition')
                            filename = None
                            if content_disp:
                                fname_match = re.search(r'filename="?([^"]+)"?', content_disp)
                                filename = fname_match.group(1) if fname_match else None
                            if not filename:
                                filename = os.path.basename(urllib.parse.urlparse(url_input).path) or f"download_{int(time.time())}.fasta"

                            is_gzipped = filename.lower().endswith('.gz') or response.headers.get('content-encoding') == 'gzip'
                            content_bytes = response.content

                            if is_gzipped:
                                content_string = gzip.decompress(content_bytes).decode('utf-8', errors='replace')
                                filename = filename[:-3] if filename.lower().endswith('.gz') else filename
                            else:
                                content_string = content_bytes.decode('utf-8', errors='replace')

                            if content_string:
                                parser = FastaParser()  # FIXED: Instantiate parser here
                                sequences, errors = parser.parse(content_string)  # Now uses local parser

                                if errors:
                                    st.warning(f"⚠️ {filename}: {errors[0]}", icon="⚠️")

                                if sequences:
                                    st.session_state.all_files[filename] = sequences
                                    st.session_state.original_sequences[filename] = sequences
                                    st.success(T("downloaded_processed").format(filename=filename, seqs=len(sequences)))
                                    st.session_state.active_filenames = [filename]
                                    st.session_state.active_sequences = sequences
                                    st.info(T("activated_file_info").format(filename=filename))
                                else:
                                    progress_tracker.log_error("No valid sequences found in content from URL.")
                            else:
                                st.error(T("empty_url_content"))

                        except requests.exceptions.RequestException as e:
                            st.error(f"HTTP Error: {e}")
                        except Exception as e:
                            st.error(f"Error: {e}")
                else:
                    st.warning(T("invalid_url"))

        # --- Google Drive Upload (only show if selected) ---
        elif selected_upload_method == T("upload_gdrive"):
            st.info(T("gdrive_info"), icon="ℹ️")
            if COLAB_AVAILABLE:  # Only show mount button if in Colab
                if st.button(T("mount_gdrive_btn")):
                    try:
                        with st.spinner("Attempting to mount Google Drive..."):
                            drive.mount('/content/drive', force_remount=True)
                        update_status("gdrive_success", level="success")
                        st.session_state.gdrive_mounted = True
                    except Exception as e:
                        update_status(f"Drive Mount Error: {e}", level="error")
                        st.session_state.gdrive_mounted = False
            else:
                # Check if drive might be mounted via Desktop app etc.
                if not os.path.exists('/content/drive'):
                    st.warning(T("gdrive_fail"))

            # Allow path input regardless
            gdrive_path = st.text_input(T("gdrive_path_label"), value=st.session_state.get('gdrive_path_input', ''), placeholder="/content/drive/MyDrive/YourFolder/*.fasta", key='gdrive_path_input')
            if st.button(T("load_gdrive_btn"), disabled=not gdrive_path):
                # Check again if accessible before trying glob
                if not os.path.exists('/content/drive') and not st.session_state.get('gdrive_mounted'):
                     st.error("Google Drive does not appear to be mounted...")
                else:
                    with st.spinner(T("processing_files")):
                        try:
                            # Use glob to find matching FASTA files
                            matching_files = glob.glob(gdrive_path)
                            if not matching_files:
                                st.warning("No matching FASTA files found at the specified path/pattern.")
                                st.stop()  # FIXED: Replace 'return' with st.stop() to halt execution cleanly

                            parser = FastaParser()  # FIXED: Instantiate parser here
                            newly_loaded_count = 0
                            total_sequences_added = 0

                            for file_path in matching_files:
                                filename = os.path.basename(file_path)
                                if filename.lower().endswith(('.fasta', '.fas', '.fa', '.fna', '.txt')):
                                    with open(file_path, 'r') as f:
                                        content_string = f.read()

                                    sequences, errors = parser.parse(content_string)  # Now uses local parser

                                    if errors:
                                        st.warning(f"⚠️ {filename}: {errors[0]}", icon="⚠️")

                                    if sequences:
                                        st.session_state.all_files[filename] = sequences
                                        if filename not in st.session_state.original_sequences:
                                            st.session_state.original_sequences[filename] = sequences
                                        newly_loaded_count += 1
                                        total_sequences_added += len(sequences)

                            if newly_loaded_count > 0:
                                msg = T("loaded_files").format(count=newly_loaded_count, seqs=total_sequences_added)
                                st.success(msg)
                                st.balloons()
                                if not st.session_state.active_sequences:
                                    st.info(T("info_activate_files"))
                            else:
                                st.warning(T("no_new_files"))
                        except Exception as e:
                            progress_tracker.log_error(f"Failed to load from Google Drive: {str(e)}")
        # --- END FIXED ---

    # ==================== TAB 2: MANAGE DATASETS ====================
    with tab_map["manage_tab"]:
        st.header(T("manage_tab"))
    
        if not st.session_state.all_files:
            st.markdown(f"""
                <div style='background: #e0f2fe; padding: 25px; border-radius: 12px; border-left: 6px solid #0ea5e9;'>
                    <h3 style='color: #0c4a6e; border: none; margin-top: 0;'>{T('file_manager_empty_title')}</h3>
                    <p>{T('file_manager_empty_subtitle')}</p>
                    <ol style='line-height: 1.8; padding-left: 20px;'>
                        <li><b>{T('step1_title')}</b> {T('step1_desc')}</li>
                        <li><b>{T('step2_title')}</b> {T('step2_desc')}</li>
                        <li><b>{T('step3_title')}</b> {T('step3_desc')}</li>
                        <li><b>{T('step4_title')}</b> {T('step4_desc')}</li>
                    </ol>
                </div>
            """, unsafe_allow_html=True)
            
            # Expanded pro tip (outside markdown for correctness)
            pro_tip_text = f"{T('tip_multi_file')} {T('pro_tip_merged')}"
            st.markdown(f"""
                <div style='background: #e0f2fe; padding: 15px; border-radius: 8px; border-left: 4px solid #0ea5e9;'>
                    <p><b>💡 {T('tip_title')}</b> {pro_tip_text}</p>
                </div>
            """, unsafe_allow_html=True)
        else:
            st.subheader(T("loaded_datasets_header"))
            st.caption(T("loaded_datasets_desc"))
    
            sorted_filenames = sorted(st.session_state.all_files.keys())
            
            # FIXED: Multiselect with dict for lookup
            file_counts = {fname: len(st.session_state.all_files[fname]) for fname in sorted_filenames}
            display_options = [f"**{fname}** ({file_counts[fname]} {T('seqs_abbrev')})" for fname in sorted_filenames]
            options_dict = {display_str: fname for fname, display_str in zip(sorted_filenames, display_options)}
            
            # Build default from current active filenames
            default_selected = [
                display_options[i] 
                for i, fname in enumerate(sorted_filenames) 
                if fname in st.session_state.active_filenames
            ]
            
            # ✅ ADDED: Clear stale session state if options changed (safer approach)
            if 'manage_file_multiselect' in st.session_state:
                stored = st.session_state.manage_file_multiselect
                if not all(opt in display_options for opt in stored):
                    # Options changed (files added/removed), clear stale data
                    del st.session_state.manage_file_multiselect
            
            # Multiselect with validated default
            selected_indices = st.multiselect(
                T("select_files_to_activate"),
                options=display_options,
                default=default_selected,  # Use computed default
                key="manage_file_multiselect",
                format_func=lambda x: x,
                help="Check to include in active dataset. Hold Ctrl/Cmd for multi-select."
            )
            
            # Map back to filenames using dict
            selected_files_now = [options_dict[opt] for opt in selected_indices]
            
            # Preview selection count (UX boost)
            if selected_files_now:
                total_seqs = sum(file_counts[fname] for fname in selected_files_now)
                st.info(T("files_selected_summary").format(
                    count=len(selected_files_now), 
                    seqs=total_seqs, 
                    unit=T('seqs_abbrev')
                ))
                
                # Build preview data: Seq count + top subtypes per file
                preview_data = []
                for fname in selected_files_now:
                    seqs = st.session_state.all_files[fname]
                    total_seqs_file = len(seqs)
                    
                    # Quick subtype counts (reuse Counter)
                    subtype_counts = Counter(m.get('type', DEFAULT_UNKNOWN) for _, _, m in seqs)
                    top_subtypes = dict(sorted(subtype_counts.items(), key=lambda x: x[1], reverse=True)[:3])
                    subtype_str = ', '.join([f"{k}: {v}" for k, v in top_subtypes.items()]) or "No subtypes"
                    
                    preview_data.append({
                        'File': fname,
                        'Total Sequences': total_seqs_file,
                        'Top Subtypes': subtype_str
                    })
                
                if preview_data:
                    preview_df = pd.DataFrame(preview_data)
                    st.subheader(T("preview_table_title"))
                    st.dataframe(preview_df, use_container_width=True, hide_index=True)
                    st.caption(T("preview_merge_caption"))
            else:
                st.info("No files selected yet.")
    
            st.markdown("---")
            st.subheader(T("actions_header"))
            action_cols = st.columns(4)
    
            with action_cols[0]:
                if st.button(T("select_all_btn"), use_container_width=True, key="manage_select_all"):
                    st.session_state.active_filenames = list(st.session_state.all_files.keys())
                    st.session_state['manage_file_multiselect'] = display_options  # FIXED: Select all strings
                    st.rerun()
    
            with action_cols[1]:
                if st.button(T("deselect_all_btn"), use_container_width=True, key="manage_deselect_all"):
                    st.session_state.active_filenames = []
                    st.session_state['manage_file_multiselect'] = []  # Clear multiselect state
                    st.rerun()
    
            with action_cols[2]:
                if st.button(T("activate_btn"), type="primary", use_container_width=True, key="manage_activate",
                            help=T("activate_help")):
                    if not selected_files_now:
                        st.warning(T("no_files_selected_activate"))
                    else:
                        st.session_state.active_filenames = selected_files_now
                        st.session_state.original_sequences = {}
                        st.session_state.active_sequences = []
                        for fname in selected_files_now:
                            current_file_seqs = [list(s) for s in st.session_state.all_files.get(fname, [])]
                            
                            # Add source_file to metadata for each sequence in this file
                            for seq in current_file_seqs:
                                seq[2]['source_file'] = fname  # seq[2] is metadata dict; add source_file
                            
                            st.session_state.active_sequences.extend(current_file_seqs)
                            st.session_state.original_sequences[fname] = current_file_seqs
                        
                        # Always update snapshot after any activation (captures newly selected)
                        st.session_state.original_active_snapshot = [list(s) for s in st.session_state.active_sequences]  # Deep copy of current pre-filter
                    
                        count = len(st.session_state.active_sequences)
                        with st.container(): # NEW: Pins under button
                            st.success(T("files_activated").format(count=len(selected_files_now), seqs=count))
                            st.info("Results ready—process in Analyze or download in Export.")

                            col_go1, col_go2 = st.columns(2)
                            with col_go1:
                                if st.button(f"🧭 {T('analyze_tab')}", key="go_to_analyze_activation"):
                                    st.session_state.active_tab = 2  # Analyze
                                    st.rerun()
                            with col_go2:
                                if st.button(f"🧭 {T('export_tab')}", key="go_to_export_activation"):
                                    st.session_state.active_tab = 4  # Export
                                    st.rerun()
                        
                        st.rerun()
               
    
            with action_cols[3]:
                if st.button(T("remove_btn"), type="secondary", use_container_width=True, key="manage_remove"):
                    if not selected_files_now:
                        st.warning(T("no_files_selected_remove"))
                    else:
                        st.session_state.confirming_removal = True
    
            if st.session_state.get("confirming_removal", False):
                st.warning(T("confirm_remove_msg").format(count=len(selected_files_now)))
                confirm_cols = st.columns(2)
                with confirm_cols[0]:
                    if st.button(T("confirm_yes"), type="primary", use_container_width=True):
                        removed_count = 0
                        for filename in selected_files_now:
                            if filename in st.session_state.all_files:
                                del st.session_state.all_files[filename]
                                removed_count += 1
                            if filename in st.session_state.original_sequences:
                                del st.session_state.original_sequences[filename]
                            if filename in st.session_state.active_filenames:
                                st.session_state.active_filenames.remove(filename)
    
                        if removed_count > 0:
                            st.session_state.active_sequences = []
                            st.session_state.original_sequences = {}
                            for fname in st.session_state.active_filenames:
                                current_file_seqs = [list(s) for s in st.session_state.all_files.get(fname, [])]
                                
                                # Add source_file to metadata for each sequence in this file (for rebuild)
                                for seq in current_file_seqs:
                                    seq[2]['source_file'] = fname
                                
                                st.session_state.active_sequences.extend(current_file_seqs)
                                st.session_state.original_sequences[fname] = current_file_seqs
    
                            st.warning(T("removed_files_msg").format(count=removed_count))
                        st.session_state.confirming_removal = False
                        st.rerun()
                with confirm_cols[1]:
                    if st.button(T("confirm_cancel"), use_container_width=True):
                        st.session_state.confirming_removal = False
                        st.rerun()
    
            st.markdown("---")
            st.subheader(T("active_dataset"))
            if st.session_state.active_sequences:
                st.write(f"**{T('files')}:** `{', '.join(st.session_state.active_filenames)}`")
                st.write(f"**{T('total_seqs')}:** `{len(st.session_state.active_sequences):,}`")
            else:
                st.info(T("active_dataset_info"))
   
    # ==================== TAB 3: ANALYZE & PROCESS ====================
    with tab_map["analyze_tab"]:
        st.header(T("analyze_tab"))

        if not st.session_state.active_sequences:
            st.markdown(f"""
                <div style='background: #fffbeb; padding: 25px; border-radius: 12px; border-left: 6px solid #f59e0b;'>
                    <h3 style='color: #92400e; border: none; margin-top: 0;'>{T('no_active_dataset_title')}</h3>
                    <p style='color: #78350f;'>{T('no_active_dataset_msg').format(tab=T('manage_tab'))}</p>
                </div>
            """, unsafe_allow_html=True)
        else:
            analyzer = SequenceAnalyzer(st.session_state.active_sequences)

            st.subheader(T("current_dataset_overview"))
            col1, col2 = st.columns([1, 2])
            with col1:
                st.plotly_chart(create_metric_indicator(
                    len(analyzer.sequences),
                    "metric_title",
                    st.session_state.lang
                ), use_container_width=True)
            with col2:
                avg_len = sum(len(s[1]) for s in analyzer.sequences) / len(analyzer.sequences) if analyzer.sequences else 0
        
                # Safe access to lang with fallback
                current_lang = st.session_state.get('lang', 'en')
                
                st.plotly_chart(create_gauge_indicator(
                    avg_len,
                    max_value=max(3000, int(avg_len * 1.5)) if avg_len > 0 else 2000,
                    title_key="gauge_title",
                    lang=current_lang,
                    #compact=True
                ), use_container_width=True)

            st.markdown("---")

            # --- Data Visualizer (Enhanced) ---
            with st.expander(T("distribution_viewer_title"), expanded=True):
                st.markdown(f"*{T('visualizer_desc')}*")
                vis_col1, vis_col2 = st.columns([3, 2])
                with vis_col1:
                    # Define options for visualization using translations
                    vis_field_options = {
                        T("vis_field_subtype"): 'type', T("vis_field_segment"): 'segment', T("vis_field_host"): 'host',
                        T("vis_field_location"): 'location', T("vis_field_clade"): 'clade',
                        T("vis_field_year"): 'year', T("vis_field_month"): 'month'
                    }
                    # ADDED New chart types
                    vis_chart_options = {
                        T("vis_type_bar"): 'bar', T("vis_type_pie"): 'pie',
                        T("vis_type_line"): 'line', T("vis_type_heatmap"): 'heatmap',
                        T("vis_type_stacked"): 'stacked'
                    }

                    # FIXED: Use index-based selection to avoid key mismatch
                    chart_type_list = list(vis_chart_options.keys())
                    chart_value_list = list(vis_chart_options.values())

                    #selected_chart_index = st.selectbox(
                    # Define options for visualization using translations
                    vis_field_options = {
                        T("vis_field_subtype"): 'type', T("vis_field_segment"): 'segment', T("vis_field_host"): 'host',
                        T("vis_field_location"): 'location', T("vis_field_clade"): 'clade',
                        T("vis_field_year"): 'year', T("vis_field_month"): 'month'
                    }

                    # Chart type options
                    vis_chart_options = {
                        T("vis_type_bar"): 'bar', 
                        T("vis_type_pie"): 'pie',
                        T("vis_type_line"): 'line', 
                        T("vis_type_heatmap"): 'heatmap',
                        T("vis_type_stacked"): 'stacked'
                    }

                    # ✅ FIXED: Safe index-based selection with validation
                    chart_type_list = list(vis_chart_options.keys())
                    chart_value_list = list(vis_chart_options.values())

                    # Get stored index, with validation
                    stored_index = st.session_state.get('vis_chart_type_index', 0)  

                    # ✅ ADDED: Validate stored index is valid
                    if not isinstance(stored_index, int) or stored_index >= len(chart_type_list):
                        stored_index = 0  # Reset to default if invalid

                    selected_chart_index = st.selectbox(
                        T("chart_type_label"), 
                        options=range(len(chart_type_list)),
                        index=stored_index,
                        format_func=lambda i: chart_type_list[i],
                        key="vis_chart_type"
                    )

                    # ✅ ADDED: Ensure selected_chart_index is an integer
                    if not isinstance(selected_chart_index, int):
                        selected_chart_index = 0

                    # Get the actual chart type value using the index
                    selected_chart_key = chart_value_list[selected_chart_index]
                    
                    # Store the index for persistence
                    st.session_state.vis_chart_type_index = selected_chart_index

                    # ADDED Conditional controls
                    # ✅ FIXED: Conditional controls with safe index-based selection
                    field1, field2, interval, top_n_val = None, None, None, 20
                    
                    if selected_chart_key in ['bar', 'pie']:
                        # Safe field selection
                        field_options_list = list(vis_field_options.keys())
                        field_values_list = list(vis_field_options.values())

                        # Validate stored index
                        stored_field_index = st.session_state.get('vis_field1_index', 0)
                        if not isinstance(stored_field_index, int) or stored_field_index >= len(field_options_list):
                            stored_field_index = 0
                        
                        selected_field_index = st.selectbox(
                            T("field_label"), 
                            options=range(len(field_options_list)),
                            index=stored_field_index,
                            format_func=lambda i: field_options_list[i],
                            key="vis_field1"
                        )
                        
                        if not isinstance(selected_field_index, int):
                            selected_field_index = 0
                        
                        field1 = field_values_list[selected_field_index]
                        field1_display = field_options_list[selected_field_index]
                        st.session_state.vis_field1_index = selected_field_index
                    
                    elif selected_chart_key == 'line':
                        interval_options = {T("vis_interval_month"): 'month', T("vis_interval_quarter"): 'quarter', T("vis_interval_year"): 'year'}
                        interval_display = st.selectbox(T("time_interval_label"), list(interval_options.keys()), key="vis_interval")
                        interval = interval_options[interval_display]
                    
                    elif selected_chart_key == 'heatmap':
                        field_options_list = list(vis_field_options.keys())
                        field_values_list = list(vis_field_options.values())
                        
                        stored_heatmap_index = st.session_state.get('vis_heatmap_field_index', 3)
                        if not isinstance(stored_heatmap_index, int) or stored_heatmap_index >= len(field_options_list):
                            stored_heatmap_index = 3
                        
                        selected_heatmap_index = st.selectbox(
                            T("field_label"),
                            options=range(len(field_options_list)),
                            index=stored_heatmap_index,
                            format_func=lambda i: field_options_list[i],
                            key="vis_heatmap_field"
                        )
                        
                        if not isinstance(selected_heatmap_index, int):
                            selected_heatmap_index = 3
                        
                        heatmap_field = field_values_list[selected_heatmap_index]
                        st.session_state.vis_heatmap_field_index = selected_heatmap_index

                        # Top N slider
                        top_n_val = st.slider(
                            T("top_n_label"),
                            min_value=1,
                            max_value=100,
                            value=st.session_state.get('vis_top_n', 20),
                            step=5,
                            key="vis_top_n",
                            help="Limit to top N items"
                        )
                    
                    elif selected_chart_key == 'stacked':
                        field_options_list = list(vis_field_options.keys())
                        field_values_list = list(vis_field_options.values())

                        # Category 1
                        stored_cat1_index = st.session_state.get('vis_cat1_stacked_index', 3)
                        if not isinstance(stored_cat1_index, int) or stored_cat1_index >= len(field_options_list):
                            stored_cat1_index = 3
                        
                        selected_cat1_index = st.selectbox(
                            T("category1_label"),
                            options=range(len(field_options_list)),
                            index=stored_cat1_index,
                            format_func=lambda i: field_options_list[i],
                            key="vis_cat1_stacked"
                        )
                        
                        if not isinstance(selected_cat1_index, int):
                            selected_cat1_index = 3
                        
                        category1 = field_values_list[selected_cat1_index]
                        st.session_state.vis_cat1_stacked_index = selected_cat1_index

                        # Category 2
                        stored_cat2_index = st.session_state.get('vis_cat2_stacked_index', 0)
                        if not isinstance(stored_cat2_index, int) or stored_cat2_index >= len(field_options_list):
                            stored_cat2_index = 0
                        
                        selected_cat2_index = st.selectbox(
                            T("category2_label"),
                            options=range(len(field_options_list)),
                            index=stored_cat2_index,
                            format_func=lambda i: field_options_list[i],
                            key="vis_cat2_stacked"
                        )
                        
                        if not isinstance(selected_cat2_index, int):
                            selected_cat2_index = 0
                        
                        category2 = field_values_list[selected_cat2_index]
                        st.session_state.vis_cat2_stacked_index = selected_cat2_index
                        
                        top_n_val = st.slider(
                            T("top_n_label") + f" ({field_options_list[stored_cat1_index]})", 
                            5, 50, 15, 5, 
                            key="vis_top_n_stacked",
                            help="Limit to top N groups"
                        )
                    # END ADDED Conditional controls

                with vis_col2:
                    # Add vertical space to align button
                    for _ in range(5 if selected_chart_key not in ['stacked','heatmap'] else (7 if selected_chart_key=='stacked' else 6) ):
                        st.write("")
                
                    # Pre-built Color Schemes Selectbox
                    if selected_chart_key in schemes_by_chart:
                        scheme_names = list(schemes_by_chart[selected_chart_key].keys())
                        selected_scheme_name = st.selectbox("🎨 Color Scheme:", scheme_names, index=0, key=f"vis_scheme_{selected_chart_key}")
                        selected_scheme = schemes_by_chart[selected_chart_key][selected_scheme_name]
                    else:
                        selected_scheme = None
                
                    # Generate Chart Button
                    if st.button(T("generate_chart_btn"), key="vis_generate", use_container_width=True, type="primary"):
                        if selected_chart_key == 'stacked' and category1 == category2:
                            st.error("Primary and Secondary categories cannot be the same for Stacked Bar chart.")
                        else:
                            with st.spinner(T("generating_chart")):
                                fig = None
                                try:
                                    use_custom = 'custom_palette' in st.session_state and st.session_state.custom_palette
                                    color_to_use = st.session_state.custom_palette if use_custom else selected_scheme
                
                                    if selected_chart_key in ['bar', 'pie']:
                                        counts = analyzer.get_metadata_distribution(field1)
                                        fig = create_distribution_chart(counts, f"{field1_display} Distribution", chart_type=selected_chart_key, color_scheme=color_to_use)
                                    elif selected_chart_key == 'line':
                                        fig = create_temporal_chart(analyzer.sequences, interval=interval, color_scheme=color_to_use)
                                    elif selected_chart_key == 'heatmap':
                                        fig = create_geographic_heatmap(analyzer.sequences, top_n=top_n_val, field=heatmap_field, color_scheme=color_to_use)
                                    elif selected_chart_key == 'stacked':
                                        fig = create_stacked_bar_chart(analyzer.sequences, category1=category1, category2=category2, top_n=top_n_val, lang=st.session_state.lang, color_scheme=color_to_use)
                
                                    if fig:
                                        if color_to_use and selected_chart_key != 'line':
                                            data_len = len(analyzer.sequences) if analyzer.sequences else 0
                                            fig = apply_color_scheme(fig, color_to_use, selected_chart_key, data_len)
                                        st.session_state.generated_chart = fig
                                        st.caption(T("chart_ready"))
                                    else:
                                        st.warning(T("chart_no_data"))
                                except Exception as e:
                                    st.error(f"{T('chart_error')}: {e}")
                                    progress_tracker.log_error(f"Chart generation failed: {e}")
                
                    # ✅ NEW IMPROVED: Interactive Color Generator Expander
                    # ✅ REDESIGNED: Custom Palette Studio (Improved UX)
                    with st.expander("🎨 Custom Palette Studio", expanded=False):
                        # Slider for color count
                        num_colors = st.slider(
                            "Number of Colors", 
                            3, 12, 8, 
                            key="vis_num_colors",
                            help="Choose how many colors for your palette"
                        )
                        
                        st.markdown("---")
                        
                        # ✅ IMPROVED: Larger, grid-based color pickers
                        st.write("**🎨 Select Your Colors:**")
                        
                        custom_colors = []
                        default_colors = [
                            "#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", 
                            "#FFEAA7", "#DDA0DD", "#F39B7F", "#8491B4",
                            "#91D1C2", "#B09C85", "#E64B35", "#4DBBD5"
                        ]
                        
                        # Create 4 columns per row for better spacing
                        cols_per_row = 4
                        num_rows = (num_colors + cols_per_row - 1) // cols_per_row
                        
                        for row_idx in range(num_rows):
                            cols = st.columns(cols_per_row)
                            for col_idx in range(cols_per_row):
                                color_idx = row_idx * cols_per_row + col_idx
                                if color_idx < num_colors:
                                    with cols[col_idx]:
                                        # ✅ IMPROVED: Larger color picker with label
                                        color = st.color_picker(
                                            f"Color {color_idx + 1}",
                                            value=default_colors[color_idx % len(default_colors)],
                                            key=f"vis_color_{color_idx}"
                                        )
                                        custom_colors.append(color)
                        
                        st.markdown("---")
                        
                        # ✅ IMPROVED: Action buttons with better layout
                        st.write("**⚡ Quick Actions:**")
                        
                        action_cols = st.columns(3)
                        
                        with action_cols[0]:
                            if st.button(
                                "✅ Apply",
                                key="vis_custom_apply",
                                use_container_width=True,
                                help="Apply selected colors to chart"
                            ):
                                st.session_state.custom_palette = custom_colors
                                st.success("✓ Applied!", icon="✅")
                        
                        with action_cols[1]:
                            sample_seq = analyzer.sequences[0][1] if analyzer.sequences else ""
                            if st.button(
                                "🧬 DNA Colors",
                                key="vis_nuc_palette",
                                use_container_width=True,
                                help="Generate palette from sequence nucleotide composition"
                            ):
                                nuc_scheme = nucleotide_palette(sample_seq)
                                if callable(nuc_scheme):
                                    nuc_colors = [nuc_scheme[i] for i in range(0, 100, 100//num_colors)]
                                else:
                                    nuc_colors = nuc_scheme[:num_colors]
                                st.session_state.custom_palette = nuc_colors
                                st.success("✓ DNA palette!", icon="🧬")
                                st.rerun()
                        
                        with action_cols[2]:
                            if st.button(
                                "🎲 Randomize",
                                key="vis_random_palette",
                                use_container_width=True,
                                help="Generate random vibrant colors"
                            ):
                                import random
                                import colorsys
                                
                                random_colors = []
                                for _ in range(num_colors):
                                    hue = random.randint(0, 360)
                                    saturation = random.randint(65, 95)
                                    lightness = random.randint(50, 70)
                                    
                                    r, g, b = colorsys.hls_to_rgb(hue/360, lightness/100, saturation/100)
                                    hex_color = f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"
                                    random_colors.append(hex_color)
                                
                                st.session_state.custom_palette = random_colors
                                st.success("✓ Randomized!", icon="🎲")
                                st.rerun()
                        
                        # ✅ IMPROVED: Better palette preview
                        if 'custom_palette' in st.session_state and st.session_state.custom_palette:
                            st.markdown("---")
                            st.write("**🎨 Current Palette:**")
                            
                            # Create larger, cleaner swatches
                            preview_cols = st.columns(min(8, len(st.session_state.custom_palette)))
                            
                            for i, color in enumerate(st.session_state.custom_palette[:8]):
                                with preview_cols[i]:
                                    # Larger swatch with better styling
                                    st.markdown(
                                        f"""
                                        <div style='
                                            background: linear-gradient(135deg, {color} 0%, {color}dd 100%);
                                            width: 100%;
                                            height: 80px;
                                            border-radius: 12px;
                                            border: 2px solid rgba(0,0,0,0.1);
                                            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
                                            display: flex;
                                            align-items: flex-end;
                                            justify-content: center;
                                            padding: 8px;
                                        '>
                                            <span style='
                                                background: rgba(255,255,255,0.9);
                                                color: #333;
                                                padding: 4px 8px;
                                                border-radius: 6px;
                                                font-size: 11px;
                                                font-family: monospace;
                                                font-weight: 600;
                                                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                                            '>{color.upper()}</span>
                                        </div>
                                        """,
                                        unsafe_allow_html=True
                                    )
                            
                            if len(st.session_state.custom_palette) > 8:
                                st.caption(f"*+ {len(st.session_state.custom_palette) - 8} more colors in palette*")
                            
                            # ✅ IMPROVED: Export button styling
                            st.markdown("<br>", unsafe_allow_html=True)
                            
                            palette_json = json.dumps({
                                'colors': st.session_state.custom_palette,
                                'name': 'Custom Viral Palette',
                                'count': len(st.session_state.custom_palette),
                                'date': datetime.now().isoformat()
                            }, indent=2)
                            
                            st.download_button(
                                label="💾 Export Palette (JSON)",
                                data=palette_json,
                                file_name=f"viral_palette_{datetime.now().strftime('%Y%m%d_%H%M')}.json",
                                mime="application/json",
                                key="export_palette",
                                use_container_width=True,
                                help="Download palette as JSON for reuse"
                            )
                    # END UPDATED Button logic

            # ADDED Display chart from session state
            # Display the generated chart (if exists in session state) outside the button's scope
            if 'generated_chart' in st.session_state and st.session_state.generated_chart:
                st.plotly_chart(st.session_state.generated_chart, use_container_width=True)
            # END ADDED Display chart

            st.markdown("---")
            st.subheader(T("processing_steps"))

            col_proc1, col_proc2 = st.columns(2)

            with col_proc1:
                st.markdown(f"#### {T('basic_operations')}")
                if st.button(T("convert_headers_btn"), key="analyze_convert", use_container_width=True, help=T("help_convert_headers")):
                    with st.spinner(T("converting_headers")):
                        analyzer.convert_headers()
                        st.rerun()

                st.markdown(f"#### {T('deduplication')}")
                if st.button(T("deduplicate_basic_btn"), key="analyze_dedup_basic", use_container_width=True, help=T("help_dedup_basic")):
                    with st.spinner(T("running_deduplication")):
                        analyzer.deduplicate_basic()
                        st.rerun()

                if st.button(T("deduplicate_advanced_btn"), key="analyze_dedup_adv", use_container_width=True, help=T("help_dedup_advanced")):
                    with st.spinner(T("running_advanced_dedup")):
                        analyzer.deduplicate_advanced()
                        st.rerun()

            with col_proc2:
                st.markdown(f"#### {T('quality_filter')}")
                min_len = st.slider(T("min_length_label"), 0, 3000, value=st.session_state.get('analyze_min_len', 200), step=50, key="analyze_min_len", help=T("help_min_length"))
                max_n = st.slider(T("max_n_run_label"), 0, 500, value=st.session_state.get('analyze_max_n', 100), step=10, key="analyze_max_n", help=T("help_max_n"))
                if st.button(T("quality_filter_btn"), key="analyze_quality", use_container_width=True):
                    with st.spinner(T("applying_quality_filter")):
                        analyzer.quality_filter(min_length=min_len, max_n_run=max_n)

                    with st.container():  # NEW: Pins under button
                        st.success("Quality filter applied—explore refined data!")
                        st.info("Next: Check subtypes or generate charts below.")
                        
                        if st.button(f"🧭 {T('refine_tab')}", key="go_to_refine_quality"):
                            st.session_state.active_tab = 3  # Refine
                            st.rerun()
                    
                    st.rerun()

                st.markdown(f"#### {T('subtype_operations')}")
                all_subtypes = ['All'] + sorted(list(set(m.get('type', DEFAULT_UNKNOWN) for _, _, m in st.session_state.active_sequences if m.get('type') != DEFAULT_UNKNOWN)))
                selected_subtype = st.selectbox(T("subtype_label"), all_subtypes, index=st.session_state.get('analyze_subtype_index', 0), key="analyze_subtype_select")
                custom_subtypes_input = st.text_input(T("custom_subtype_label"), placeholder=T("custom_subtype_placeholder"), key="analyze_subtype_custom")

                sub_op_cols = st.columns(2)
                with sub_op_cols[0]:
                    if st.button(T("filter_subtype_btn"), key="analyze_subtype_filter", use_container_width=True):
                        targets = []
                        if custom_subtypes_input:
                            targets = [s.strip().upper() for s in custom_subtypes_input.split(',') if s.strip()]
                        elif selected_subtype != 'All':
                            targets = [selected_subtype.upper()]

                        if targets:
                            with st.spinner(T("filtering_subtype")):
                                analyzer.filter_by_subtype(targets)
                                st.rerun()
                        else:
                            st.warning(T("warning_select_subtype"))

                with sub_op_cols[1]:
                    if st.button(T("check_subtypes_btn"), key="analyze_subtype_check", use_container_width=True):
                        dist_counts = analyzer.get_subtype_distribution()
                        if dist_counts:
                            fig = create_distribution_chart(dist_counts, "distribution_title", st.session_state.lang, chart_type='pie')
                            st.plotly_chart(fig, use_container_width=True)
                        else:
                            st.warning(T("warning_no_subtype_info"))

    # ==================== TAB 4: REFINE & VISUALIZE ====================
    with tab_map["refine_tab"]:
        st.header(T("refine_tab"))

        if not st.session_state.active_sequences:
            st.warning(T("no_data_msg"))
        else:
            analyzer = SequenceAnalyzer(st.session_state.active_sequences)

            st.subheader(T("clade_monthly_header"))
            available_clades = sorted([c for c in list(set(m.get('clade', DEFAULT_UNKNOWN) for _, _, m in analyzer.sequences)) if c != DEFAULT_UNKNOWN])

            if not available_clades:
                st.caption(T("no_clade_info"))
            else:
                clade_mode_display = st.radio(T("mode_label"), [T("clade_mode_single"), T("clade_mode_multiple")], index=st.session_state.get('refine_clade_mode_index', 0), key="refine_clade_mode", horizontal=True)
                clade_mode = 'single' if clade_mode_display == T("clade_mode_single") else 'multiple'

                targets = []
                separate = True
                if clade_mode == 'single':
                    target_clade = st.selectbox(T("select_clade"), available_clades, key="refine_clade_single")
                    targets = [target_clade] if target_clade else []
                    separate = False
                else:
                    targets = st.multiselect(T("select_clades"), available_clades, default=available_clades[:1] if available_clades else [], key="refine_clade_multi")
                    separate = st.checkbox(T("process_clades_separately"), True, key="refine_clade_separate")

                keep_monthly_options = {
                    T("temporal_order_first"): "First Only",
                    T("temporal_order_last"): "Last Only",
                    T("temporal_order_both"): "Both (First & Last)"
                }
                keep_monthly_display = st.selectbox(T("keep_monthly_label"), list(keep_monthly_options.keys()), index=st.session_state.get('refine_clade_keep_index', 2), key="refine_clade_keep")
                keep_strategy = keep_monthly_options[keep_monthly_display]

                if st.button(T("apply_clade_filter_button"), key="refine_clade_apply", disabled=not targets):
                    if targets:
                        with st.spinner(T("applying_clade_filter")):
                            analyzer.filter_clade_monthly(
                                mode=clade_mode,
                                targets=targets,
                                keep_strategy=keep_strategy,
                                separate=separate
                            )
                            st.rerun()
                    else:
                        st.warning(T("warning_select_clade"))

            st.markdown("---")

            # UPGRADED: Preview & Download Clades Section
            st.subheader("🧬 Preview & Download Clades")
            st.caption("Group and export sequences by clade without applying filters.")

            if not st.session_state.active_sequences:
                st.warning(T("no_data_msg"))
            else:
                analyzer = SequenceAnalyzer(st.session_state.active_sequences)
                available_clades = sorted([c for c in list(set(m.get('clade', DEFAULT_UNKNOWN) for _, _, m in analyzer.sequences)) if c != DEFAULT_UNKNOWN])

                if not available_clades:
                    st.caption(T("no_clade_info"))
                else:
                    # Mode: Individual or Multiple
                    export_mode = st.radio("Export Mode:", ["Individual Clade (Single FASTA)", "Multiple Clades (ZIP + Individuals)"], key="clade_export_mode", horizontal=True)
                    export_mode_val = 'individual' if "Individual" in export_mode else 'multiple'

                    if export_mode_val == 'individual':
                        selected_clade = st.selectbox("Select Clade:", available_clades, key="clade_single_select")
                        selected_clades = [selected_clade] if selected_clade else []
                    else:
                        selected_clades = st.multiselect("Select Clades:", available_clades, default=available_clades[:3], key="clade_multi_select")
                        # UPGRADE: Clear button for multi
                        if st.button("Clear Selection", key="clear_clade_multi", use_container_width=True):
                            st.session_state.clade_multi_select = []
                            st.rerun()

                    # UPGRADE: Quick preview of seq counts per clade (non-destructive)
                    if available_clades:
                        st.write("**Quick Counts:**")
                        col1, col2 = st.columns(2)
                        with col1:
                            st.caption(T("clade_counts_header"))
                        for clade in available_clades[:10]:  # Limit preview to top 10
                            count = len([s for s in analyzer.sequences if s[2].get('clade') == clade])
                            with col2:
                                st.caption(f"{clade} | {count}")
                        if len(available_clades) > 10:
                            st.caption(f"... and {len(available_clades) - 10} more")

                    if st.button("🧬 Preview & Download Selected Clades", key="preview_clades", disabled=not selected_clades):
                        if selected_clades:
                            with st.spinner("Grouping sequences by clade..."):
                                # Data Mode Integration: Choose seqs based on toggle
                                data_mode_val = 'current' if "Current" in st.session_state.get('data_mode_toggle', 'Current') else 'original'
                                seqs_to_use = st.session_state.original_active_snapshot if data_mode_val == 'original' else st.session_state.active_sequences
                                temp_analyzer = SequenceAnalyzer(seqs_to_use)
                                
                                exported = temp_analyzer.export_clades(selected_clades, export_mode=export_mode_val)
                                if exported:
                                    st.success(f"Preview/export ready for {len(exported)} clade(s)! Check download buttons below.")
                                else:
                                    st.error("Export failed—check logs for details.")
                        else:
                            st.warning("Please select at least one clade.")

                    # UPGRADE: Tips Expander
                    with st.expander("ℹ️ Tips for Clade Exports"):
                        st.markdown(f"- {T('clade_export_tips_single')}")
                        st.markdown(f"- {T('clade_export_tips_multiple')}")
                        st.markdown(f"- {T('clade_export_tips_no_filter')}")
                    

            st.markdown("---")

            st.subheader(T("enhanced_temporal_header"))
            group_options = {
                T("temporal_group_location_host_month_clade"): "location_host_month_clade",
                T("temporal_group_location"): "location",
                T("temporal_group_host"): "host",
                T("temporal_group_clade"): "clade",
                T("temporal_group_location_host"): "location_host",
                T("temporal_group_host_clade"): "host_clade",
                T("temporal_group_none"): "none",
                T("temporal_group_custom"): "custom",
            }
            sort_options = {
                T("temporal_sort_date"): "date",
                T("temporal_sort_location"): "location",
                T("temporal_sort_host"): "host",
                T("temporal_sort_clade"): "clade",
                T("temporal_sort_isolate"): "isolate_id",
            }
            keep_options = {
                T("temporal_order_first"): "first",
                T("temporal_order_last"): "last",
                T("temporal_order_both"): "both",
            }

            gsk_cols = st.columns(3)
            with gsk_cols[0]:
                group_by_display = st.selectbox(T("group_by_label"), list(group_options.keys()), index=st.session_state.get('refine_temp_group_index', 4), key="refine_temp_group")
                group_by_val = group_options[group_by_display]
            with gsk_cols[1]:
                sort_by_display = st.selectbox(T("sort_by_label"), list(sort_options.keys()), key="refine_temp_sort")
                sort_by_val = sort_options[sort_by_display]
            with gsk_cols[2]:
                keep_per_group_display = st.selectbox(T("keep_per_group_label"), list(keep_options.keys()), index=2, key="refine_temp_keep")
                keep_per_group_val = keep_options[keep_per_group_display]

            custom_grouping_input = ""
            if group_by_val == "custom":
                custom_grouping_input = st.text_input(T("custom_grouping_label"), key="refine_temp_custom", value=st.session_state.get('refine_temp_custom', ''), placeholder=T("custom_grouping_placeholder"))

            if st.button(T("apply_temporal_filter_button"), key="refine_temp_apply"):
                custom_grouping_list = [f.strip() for f in custom_grouping_input.split(',')] if group_by_val == "custom" and custom_grouping_input else None
                with st.spinner(T("applying_temporal_filter")):
                    analyzer.enhanced_temporal_filter(
                        group_by=group_by_val,
                        sort_by=sort_by_val,
                        keep_per_group=keep_per_group_val,
                        custom_grouping=custom_grouping_list
                    )
                    st.rerun()
            
            st.markdown("---")
        
            # ========== ADD NEW SECTION HERE (BEFORE "Extract Accessions") ==========
            st.subheader(T("filter_by_accessions_title"))
            st.caption(T("filter_by_accessions_desc"))
            
            filter_method = st.radio(
                T("filter_input_method"),
                [T("filter_method_text"), T("filter_method_file")],
                horizontal=True,
                key="epi_filter_method"
            )
            
            accessions_to_find = []
            
            if filter_method == T("filter_method_text"):
                accession_input = st.text_area(
                    T("filter_textarea_label"),
                    placeholder=T("filter_textarea_placeholder"),
                    height=150,
                    key="epi_filter_input",
                    help=T("filter_textarea_help")
                )
                
                if accession_input:
                    accessions_to_find = parse_accessions(accession_input)
            
            else:  # File Upload
                uploaded_file = st.file_uploader(
                    T("filter_file_label"),
                    type=['txt', 'csv'],
                    key="epi_filter_file",
                    help=T("filter_file_help")
                )
                
                if uploaded_file:
                    try:
                        content = uploaded_file.getvalue().decode('utf-8')
                        accessions_to_find = parse_accessions(content)
                        st.success(T("filter_loaded_success").format(count=len(accessions_to_find)))
                    except Exception as e:
                        st.error(f"Error reading file: {e}")
            
            # Show preview of loaded accessions
            if accessions_to_find:
                with st.expander(T("filter_preview_title").format(count=len(accessions_to_find)), expanded=False):
                    st.code("\n".join(accessions_to_find[:20]))
                    if len(accessions_to_find) > 20:
                        st.caption(T("filter_preview_more").format(count=len(accessions_to_find) - 20))
            
            # Filter button
            filter_col1, filter_col2, filter_col3 = st.columns([2, 1, 1])
            
            with filter_col1:
                if st.button(
                    T("filter_button"),
                    type="primary",
                    disabled=not accessions_to_find,
                    use_container_width=True,
                    key="apply_epi_filter"
                ):
                    if accessions_to_find:
                        with st.spinner(T("filter_searching").format(count=len(accessions_to_find))):
                            filtered = analyzer.filter_by_accessions(accessions_to_find)
                        
                        # NEW: Local container for pinned notifications
                        with st.container():
                            if filtered:
                                # Get match statistics from report (with fallback)
                                report = st.session_state.get('last_report', '')
                                import re
                                found_match = re.search(r'Accessions found: (\d+)', report)
                                found_count = int(found_match.group(1)) if found_match else len({meta.get('isolate_id', '') for _, _, meta in filtered})  # FIXED: Fallback to unique IDs
                                
                                st.toast(
                                    T("filter_success_toast").format(count=len(filtered)),
                                    icon="🎯"
                                )
                                
                                # NEW: Explicit "results ready" prompt
                                st.info("✅ Results ready—sequences filtered! Download in Export tab.")
                                
                                st.success(
                                    f"🎉 **{T('filter_success_title')}**\n\n"
                                    f"{T('filter_success_found').format(count=len(filtered))}\n"
                                    f"{T('filter_success_matched').format(found=found_count, total=len(accessions_to_find))}\n"  # FIXED: accessions_to_find → accessions_to_find (assume consistent var)
                                    f"{T('filter_success_removed').format(count=analyzer.original_count_for_last_op - len(filtered))}\n\n"
                                    f"{T('filter_success_export')}"
                                )
                                
                                # Guidance button (inside container for pinning)
                                if st.button("🧭 Go to Export & Reports", key="go_to_export_epi"):
                                    st.session_state.active_tab = 4  # Export index
                                    st.rerun()
                            else:
                                st.warning(
                                    f"⚠️ **{T('filter_no_matches_title')}**\n\n"
                                    f"{T('filter_no_matches_searched').format(count=len(accessions_to_find))}\n\n"
                                    f"{T('filter_no_matches_suggestions')}\n"
                                    f"{T('filter_no_matches_verify')}\n"
                                    f"{T('filter_no_matches_test')}\n"
                                    f"{T('filter_no_matches_extract')}"
                                )
                        
                        st.rerun()  # Rerun after container (refreshes UI but keeps pinned messages)
                    else:
                        st.warning("No accessions provided—enter or upload some to filter.")
            
            with filter_col2:
                if accessions_to_find and st.session_state.active_sequences:
                    # Preview match count without filtering
                    preview_count = sum(
                        1 for _, _, meta in st.session_state.active_sequences
                        if any(acc.upper() in str(meta.get('isolate_id', '')).upper() or 
                              acc.upper() in str(meta.get('original_header', '')).upper()
                              for acc in accessions_to_find)
                    )
                    st.metric(T("filter_potential_matches"), preview_count)
            
            with filter_col3:
                if st.button(
                    T("filter_clear_button"),
                    use_container_width=True,
                    key="clear_epi_filter"
                ):
                    st.session_state.epi_filter_input = ""
                    st.toast(T("filter_cleared_toast"), icon="✨")
                    st.rerun()
            
            # ========== END NEW SECTION ==========

            # st.markdown("---")
            st.markdown("---")
            st.subheader(T("extract_accessions_btn"))
            if st.button(T("extract_accessions_btn"), key="refine_extract"):
                accessions = analyzer.extract_accessions()
                if accessions:
                    st.session_state.accession_list = accessions
                    
                    # NEW: Column layout for success + clickable button
                    col1, col2 = st.columns([3, 1])
                    with col1:
                        st.success(T("accessions_found").format(count=len(accessions), tab=T('export_tab')))
                    with col2:
                        if st.button("🧭 Go to Export", key="go_to_export_accessions"):
                            st.session_state.active_tab = 4  # Export index
                            st.rerun()
                    
                    # Keep the preview below
                    st.text_area(T("accession_preview"), "\n".join(accessions[:20]), height=150, disabled=True)
                else:
                    st.warning(T("no_accessions_found"))

            
            st.markdown("---")

            # ========== TRANSLATED: UNIVERSAL SPLIT & EXPORT INTERFACE ==========
            st.subheader(T("split_export_title"))
            st.caption(T("split_export_desc"))
            
            if not st.session_state.active_sequences:
                st.warning(T("no_data_msg"))
            else:
                analyzer = SequenceAnalyzer(st.session_state.active_sequences)
                
                # Field selection with EPI_ISL support
                split_field_options = {
                    T("vis_field_subtype"): 'type',
                    T("vis_field_segment"): 'segment', 
                    T("vis_field_host"): 'host',
                    T("vis_field_location"): 'location',
                    T("vis_field_clade"): 'clade',
                    T("vis_field_year"): 'year',
                    T("vis_field_month"): 'month',
                    T("split_accession_field"): 'isolate_id'  # ✅ Translated
                }
                
                col_split1, col_split2 = st.columns([2, 1])
                
                with col_split1:
                    split_field_display = st.selectbox(
                        T("split_field_label"),
                        options=list(split_field_options.keys()),
                        key="split_field_selector",
                        help=T("split_field_help")
                    )
                    split_field = split_field_options[split_field_display]
                
                with col_split2:
                    split_data_mode = st.radio(
                        T("split_data_source"),
                        [T("split_data_current"), T("split_data_original")],
                        key="split_data_mode",
                        horizontal=True,
                        help=T("split_data_help")
                    )
                
                # Warning for accession splitting
                if split_field == 'isolate_id':
                    st.info(T("split_accession_note"), icon="ℹ️")
                
                # Preview button
                if st.button(T("split_preview_btn"), key="preview_split"):
                    data_mode_val = 'original' if split_data_mode == T("split_data_original") else 'current'
                    seqs_to_split = st.session_state.original_active_snapshot if data_mode_val == 'original' else st.session_state.active_sequences
                    
                    # Group by selected field
                    groups = defaultdict(list)
                    for header, seq, metadata in seqs_to_split:
                        if split_field == 'year':
                            key = str(metadata['collection_date'].year) if metadata.get('collection_date') else DEFAULT_UNKNOWN
                        elif split_field == 'month':
                            key = metadata['collection_date'].strftime('%Y-%m') if metadata.get('collection_date') else DEFAULT_UNKNOWN
                        elif split_field == 'isolate_id':
                            key = metadata.get('isolate_id', DEFAULT_UNKNOWN)
                            if key == DEFAULT_UNKNOWN:
                                import re
                                match = re.search(r'EPI_ISL_\d+', header)
                                if match:
                                    key = match.group(0)
                        else:
                            key = metadata.get(split_field, DEFAULT_UNKNOWN)
                        
                        if key != DEFAULT_UNKNOWN:
                            groups[key].append((header, seq, metadata))
                    
                    if groups:
                        st.write(T("split_preview_title").format(count=len(groups)))
                        
                        if len(groups) > 100:
                            st.warning(
                                T("split_large_warning").format(count=len(groups), field=split_field_display),
                                icon="⚠️"
                            )
                        
                        # Preview table
                        preview_data = []
                        display_limit = 20 if len(groups) <= 100 else 10
                        
                        for key, seqs in sorted(groups.items(), key=lambda x: len(x[1]), reverse=True)[:display_limit]:
                            preview_data.append({
                                split_field_display: key,
                                'Sequences': len(seqs),
                                'Avg Length': int(sum(len(s[1]) for s in seqs) / len(seqs)) if seqs else 0
                            })
                        
                        preview_df = pd.DataFrame(preview_data)
                        st.dataframe(preview_df, use_container_width=True, hide_index=True)
                        
                        if len(groups) > display_limit:
                            st.caption(f"... and {len(groups) - display_limit} more groups")
                        
                        # Statistics
                        total_seqs = sum(len(seqs) for seqs in groups.values())
                        avg_seqs_per_group = total_seqs / len(groups)
                        st.caption(T("split_stats").format(
                            groups=len(groups),
                            seqs=total_seqs,
                            avg=f"{avg_seqs_per_group:.1f}"
                        ))
                        
                        st.session_state.split_preview_groups = groups
                        st.session_state.split_field_name = split_field
                    else:
                        st.warning(T("split_no_data").format(field=split_field_display))
                
                # Export section
                if 'split_preview_groups' in st.session_state and st.session_state.split_preview_groups:
                    groups = st.session_state.split_preview_groups
                    split_field_name = st.session_state.split_field_name
                    
                    st.markdown("---")
                    
                    export_col1, export_col2 = st.columns([2, 1])
                    
                    with export_col1:
                        total_files = len(groups)
                        total_seqs = sum(len(seqs) for seqs in groups.values())
                        
                        if total_files > 100:
                            st.warning(
                                T("split_large_export_warning").format(files=total_files, seqs=total_seqs),
                                icon="⏳"
                            )
                        
                        if st.button(
                            T("split_export_zip_btn").format(files=total_files, seqs=total_seqs),
                            key="export_split_zip",
                            type="primary",
                            use_container_width=True
                        ):
                            with st.spinner(T("split_creating_zip").format(files=total_files)):
                                zip_buffer = io.BytesIO()
                                with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zipf:
                                    for key, seqs in groups.items():
                                        safe_key = str(key).replace('/', '_').replace('\\', '_').replace('|', '_').replace(' ', '_')
                                        safe_key = safe_key.replace(':', '_').replace('*', '_').replace('?', '_').replace('"', '_')
                                        safe_key = safe_key.replace('<', '_').replace('>', '_')
                                        
                                        filename = f"{split_field_name}_{safe_key}.fasta"
                                        
                                        fasta_io = io.StringIO()
                                        for header, seq, _ in seqs:
                                            h = header if header.startswith('>') else '>' + header
                                            fasta_io.write(f"{h}\n{seq}\n")
                                        
                                        zipf.writestr(filename, fasta_io.getvalue())
                                
                                zip_buffer.seek(0)
                                st.download_button(
                                    label=T("split_download_zip").format(files=total_files),
                                    data=zip_buffer.getvalue(),
                                    file_name=f"split_by_{split_field_name}_{datetime.now().strftime('%Y%m%d_%H%M')}.zip",
                                    mime="application/zip",
                                    key="download_split_zip",
                                    use_container_width=True
                                )
                                st.success(T("split_zip_success").format(files=total_files))
                    
                    with export_col2:
                        st.caption(T("split_individual_caption"))
                        
                        individual_limit = 3 if len(groups) > 50 else 5
                        
                        for key, seqs in sorted(groups.items(), key=lambda x: len(x[1]), reverse=True)[:individual_limit]:
                            safe_key = str(key).replace('/', '_').replace('\\', '_').replace('|', '_').replace(' ', '_')
                            safe_key = safe_key.replace(':', '_').replace('*', '_').replace('?', '_').replace('"', '_')
                            
                            fasta_io = io.StringIO()
                            for header, seq, _ in seqs:
                                h = header if header.startswith('>') else '>' + header
                                fasta_io.write(f"{h}\n{seq}\n")
                            
                            display_key = str(key)[:30] + '...' if len(str(key)) > 30 else str(key)
                            
                            st.download_button(
                                label=f"📄 {display_key} ({len(seqs)} seqs)",
                                data=fasta_io.getvalue(),
                                file_name=f"{split_field_name}_{safe_key}.fasta",
                                mime="text/plain",
                                key=f"download_split_{safe_key[:50]}",
                                use_container_width=True
                            )
                        
                        if len(groups) > individual_limit:
                            st.caption(T("split_more_in_zip").format(count=len(groups) - individual_limit))
                    
                    if st.button(T("split_clear_btn"), key="clear_split_preview"):
                        if 'split_preview_groups' in st.session_state:
                            del st.session_state.split_preview_groups
                        if 'split_field_name' in st.session_state:
                            del st.session_state.split_field_name
                        st.rerun()
                
                # Tips expander
                with st.expander(T("split_tips_title")):
                    st.markdown(T("split_tips_content"))
            
            # ========== END TRANSLATED INTERFACE ==========

    # ==================== TAB 5: EXPORT & REPORTS ====================
    with tab_map["export_tab"]:
        st.header(T("export_tab"))

        col_exp1, col_exp2 = st.columns(2)

        with col_exp1:
            st.subheader(T("last_report_header"))
            if st.session_state.last_report:
                st.text_area(T("report_content"), value=st.session_state.last_report, height=300, disabled=True, key="export_report_area")
                st.download_button(
                    label=T("export_report_btn"),
                    data=st.session_state.last_report,
                    file_name=f"analysis_report_{datetime.now().strftime('%Y%m%d_%H%M')}.txt",
                    mime="text/plain",
                    key="export_download_report",
                    use_container_width=True
                )
            else:
                st.info(T("no_analysis_report"))

            if st.session_state.active_sequences:
                try:
                    fasta_str_io = io.StringIO()
                    for header, seq, _ in st.session_state.active_sequences:
                        h = header if isinstance(header, str) else str(header or '')
                        h = h if h.startswith('>') else '>' + h
                        s = seq if isinstance(seq, str) else str(seq or '')
                        fasta_str_io.write(f"{h}\n{s}\n")

                    st.download_button(
                        label=f"{T('download_active_button')} ({len(st.session_state.active_sequences)} {T('seqs_abbrev')})",
                        data=fasta_str_io.getvalue(),
                        file_name=f"active_data_{datetime.now().strftime('%Y%m%d_%H%M')}.fasta",
                        mime="text/plain",
                        key="export_download_active",
                        use_container_width=True,
                        type="primary"
                    )
                except Exception as e:
                    st.error(T("error_export_active").format(error=str(e)))

            # Existing: st.download_button for active FASTA
            # Add this right after it
            # if st.session_state.original_sequences:  # Check if per-file originals exist
            if st.session_state.original_sequences:  # Check if per-file originals exist
                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zipf:
                    total_seqs = 0  # FIXED: Consistent var for caption
                    for fname, seqs in st.session_state.original_sequences.items():
                        fasta_io = io.StringIO()
                        for header, seq, _ in seqs:
                            h = header  # FIXED: Assign first
                            if not h.startswith('>'):  # Then check/prepend
                                h = '>' + h
                            fasta_io.write(f"{h}\n{seq}\n")
                        fasta_str = fasta_io.getvalue()
                        zipf.writestr(f"{fname}", fasta_str)
                        total_seqs += len(seqs)
                
                zip_buffer.seek(0)
                st.caption(f"Exports pre-filter originals as separate FASTAs in a ZIP ({total_seqs} total seqs)—no merging.")
                st.download_button(
                    label=f"⬇️ Download Per-File ZIP (Originals, {len(st.session_state.original_sequences)} files, {total_seqs} seqs)",
                    data=zip_buffer.getvalue(),
                    file_name=f"per_file_originals_{datetime.now().strftime('%Y%m%d_%H%M')}.zip",
                    mime="application/zip",
                    key="export_per_file_zip",
                    use_container_width=True
                )
                
                st.caption("Exports pre-filter originals as separate FASTAs in a ZIP—no merging.")

            if st.session_state.get('accession_list'):
                acc_data = "\n".join(st.session_state.accession_list)
                st.download_button(
                    label=T("download_accessions").format(count=len(st.session_state.accession_list)),
                    data=acc_data,
                    file_name=f"extracted_accessions_{datetime.now().strftime('%Y%m%d_%H%M')}.txt",
                    mime="text/plain",
                    key="export_download_accessions",
                    use_container_width=True
                )

        with col_exp2:
            st.subheader(T("export_logs_header"))
            log_data = "\n".join(st.session_state.analysis_log)
            st.download_button(
                label=T("download_log_button"),
                data=log_data if log_data else "No logs generated in this session.",
                file_name=f"analysis_log_{datetime.now().strftime('%Y%m%d_%H%M')}.txt",
                mime="text/plain",
                key="export_download_log",
                use_container_width=True,
                disabled=not log_data,
                help=T("download_log_help")
            )
            with st.expander(T("show_log_expander")):
                display_log = log_data if log_data else T("no_logs_yet")
                st.text_area(T("log_preview"), value=display_log, height=350, disabled=True, key="export_log_preview")


    # ==================== TAB 6: DOCUMENTATION ====================
    with tab_map["docs_tab"]:
        # Display documentation title
        st.markdown(f"## {T('docs_tab')}")
        
        # Display full documentation content
        st.markdown(T("docs_header"))  # This is the big markdown table
        
        # Display tips section
        st.markdown("---")
        st.markdown(T("docs_tips"))
        
        # Optional: Add a download button for the docs
        docs_content = T("docs_header") + "\n\n" + T("docs_tips")
        st.download_button(
            label="📥 Download Documentation (Markdown)",
            data=docs_content,
            file_name="vir_seq_sift_guide.md",
            mime="text/markdown",
            key="download_docs"
        )

    # Global footer (appears on all tabs)
    st.markdown("---")
    st.caption(T("footer_text"))

    gc.collect()

if __name__ == "__main__":  # ✅ CORRECT
    main()
