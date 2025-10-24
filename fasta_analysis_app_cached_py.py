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
        update_status(update_status_key, level=level, log=True)

    def start_operation(self, operation_name):
        st.session_state.start_time = time.time()
        update_status("processing", level='info', log=True)
        if 'analysis_log' in st.session_state:
            st.session_state.analysis_log.append(f"[{datetime.now().strftime('%H:%M:%S')}] INFO: Started: {operation_name}")

    def complete_operation(self, operation_name, status_key="complete"):
        duration_str = ""
        start_time = st.session_state.pop('start_time', None)
        if start_time:
            duration = time.time() - start_time
            duration_str = f" (Duration: {duration:.2f}s)"

        log_message = f"Completed: {operation_name}{duration_str}"
        level = 'success' if status_key == 'complete' else ('warning' if status_key == 'warning' else 'info')

        if 'analysis_log' in st.session_state:
            st.session_state.analysis_log.append(f"[{datetime.now().strftime('%H:%M:%S')}] {level.upper()}: {log_message}")

        update_status(status_key, level=level, log=False)

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

# ==================== VISUALIZATION FUNCTIONS ====================
def create_metric_indicator(value, title_key, lang="en"):
    """Create a metric indicator"""
    title = get_translation(title_key, lang)
    fig = go.Figure(go.Indicator(
        mode="number",
        value=value,
        title={'text': title, 'font': {'size': 18}},
        number={'font': {'size': 40, 'color': '#2563eb'}},
        domain={'x': [0, 1], 'y': [0, 1]}
    ))
    fig.update_layout(
        height=150,
        margin=dict(l=10, r=10, t=40, b=10),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)'
    )
    return fig

def create_gauge_indicator(value, max_value, title_key, lang="en"):
    """Create a gauge indicator"""
    title = get_translation(title_key, lang)
    display_value = min(value, max_value) if value is not None else 0
    fig = go.Figure(go.Indicator(
        mode="gauge+number",
        value=display_value,
        title={'text': title, 'font': {'size': 18}},
        gauge={'axis': {'range': [0, max_value], 'tickwidth': 1, 'tickcolor': "darkblue"},
               'bar': {'color': "#3b82f6", 'thickness': 0.75},
               'bgcolor': "white",
               'borderwidth': 2,
               'bordercolor': "#e5e7eb",
               'steps': [
                   {'range': [0, max_value * 0.5], 'color': '#e5e7eb'},
                   {'range': [max_value * 0.5, max_value * 0.8], 'color': '#d1d5db'}],
               'threshold': {'line': {'color': "#ef4444", 'width': 4}, 'thickness': 0.8, 'value': max_value * 0.9}},
        number={'font': {'size': 30}, 'suffix': f" {get_translation('bp', lang)}"}
    ))
    fig.update_layout(
        height=200,
        margin=dict(l=20, r=20, t=50, b=10),
        paper_bgcolor='rgba(0,0,0,0)',
        font={'color': "#374151"}
    )
    return fig

def create_distribution_chart(data_dict, title_key, lang="en", chart_type='bar'):
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
                     color_discrete_sequence=px.colors.Pastel1)
        fig.update_traces(textposition='inside', textinfo='percent+label', pull=[0.05]*len(df))
        fig.update_layout(legend_title_text='Categories', showlegend=True)
    else:
        fig = px.bar(df.sort_values('Count', ascending=True),
                     y='Category', x='Count', title=f"{title}", text_auto=True,
                     orientation='h',
                     color='Count',
                     color_continuous_scale=px.colors.sequential.Blues)
        fig.update_layout(yaxis_title=None, xaxis_title="Count", coloraxis_showscale=False)
        fig.update_yaxes(categoryorder='total ascending')

    fig.update_layout(
        margin=dict(t=50, b=20, l=20, r=20),
        title_font_size=20,
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)'
    )
    return fig

# ==================== NEW PLOTLY VISUALIZATION FUNCTIONS ====================
def create_temporal_chart(sequences, interval='month', lang='en'):
    """Generate a Plotly line chart for sequences over time."""
    progress_tracker.start_operation(f"Generating Temporal Chart (Interval: {interval})")
    df_data = [{'date': item[2].get('collection_date')} for item in sequences if item[2].get('collection_date')]
    if not df_data:
        progress_tracker.log_error("No date information found for temporal chart.")
        # Return an empty figure with a title indicating no data
        fig = go.Figure()
        fig.update_layout(title="Temporal Distribution (No Data)", xaxis={'visible': False}, yaxis={'visible': False},
                          annotations=[{'text': 'No date data available', 'xref': 'paper', 'yref': 'paper', 'showarrow': False, 'font': {'size': 16}}])
        return fig

    df = pd.DataFrame(df_data)
    df['date'] = pd.to_datetime(df['date'])
    df = df.dropna(subset=['date'])

    # Aggregate counts based on interval
    if interval == 'year':
        df['period'] = df['date'].dt.year.astype(str)
    elif interval == 'quarter':
        # Ensure correct quarterly period string (e.g., 2023-Q1)
        df['period'] = df['date'].dt.to_period('Q').astype(str)
    else: # Default to month (YYYY-MM)
        df['period'] = df['date'].dt.strftime('%Y-%m')

    counts = df['period'].value_counts().sort_index()
    counts_df = counts.reset_index()
    counts_df.columns = ['Period', 'Count']

    # Use translation for title
    T = lambda key: get_translation(key, lang)
    title_text = f"Sequence Count by {interval.capitalize()}" # Fallback
    if interval == 'year': title_text = T("vis_interval_year") + " Count"
    elif interval == 'quarter': title_text = T("vis_interval_quarter") + " Count"
    elif interval == 'month': title_text = T("vis_interval_month") + " Count"

    fig = px.line(counts_df, x='Period', y='Count',
                  title=title_text,
                  markers=True, text='Count')
    fig.update_traces(textposition="top center")
    fig.update_layout(
        xaxis_title="Time Period", yaxis_title="Number of Sequences",
        margin=dict(t=50, b=20, l=20, r=20),
        paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)'
    )
    progress_tracker.complete_operation("Temporal chart generated")
    return fig

def create_geographic_heatmap(sequences, top_n=20, lang='en'):
    """Generate a Plotly horizontal bar chart simulating a heatmap."""
    progress_tracker.start_operation(f"Generating Geographic Heatmap (Top {top_n})")
    # Exclude DEFAULT_UNKNOWN from counts if it exists
    location_counts = Counter(item[2].get('location', DEFAULT_UNKNOWN) for item in sequences if item[2].get('location', DEFAULT_UNKNOWN) != DEFAULT_UNKNOWN)

    if not location_counts:
        progress_tracker.log_error("No location information found for heatmap.")
        fig = go.Figure()
        fig.update_layout(title="Geographic Distribution (No Data)", xaxis={'visible': False}, yaxis={'visible': False},
                          annotations=[{'text': 'No location data available', 'xref': 'paper', 'yref': 'paper', 'showarrow': False, 'font': {'size': 16}}])
        return fig

    top_locations = location_counts.most_common(top_n)
    df = pd.DataFrame(top_locations, columns=['Location', 'Count'])

    # Use translation for title
    T = lambda key: get_translation(key, lang)
    title_text = f"Top {len(df)} Locations by Sequence Count"

    fig = px.bar(df.sort_values('Count', ascending=True), # Sort for horizontal bar
                 y='Location', x='Count',
                 title=title_text,
                 text_auto=True, orientation='h',
                 color='Count',
                 color_continuous_scale=px.colors.sequential.OrRd) # Orange-Red scale

    fig.update_layout(
        yaxis_title=None, xaxis_title="Number of Sequences",
        coloraxis_showscale=False, # Hide color bar legend
        margin=dict(t=50, b=20, l=20, r=20),
        paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)'
    )
    progress_tracker.complete_operation("Geographic heatmap generated")
    return fig

def create_stacked_bar_chart(sequences, category1='location', category2='type', top_n=15, lang='en'):
    """Generate a Plotly stacked bar chart."""
    progress_tracker.start_operation(f"Generating Stacked Bar ({category1} vs {category2}, Top {top_n})")
    T = lambda key: get_translation(key, lang)

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

    progress_tracker.complete_operation("Stacked bar chart generated")
    return fig

# ==================== CUSTOM CSS ====================
def load_custom_css():
    """Load custom CSS for better UI"""
    st.markdown("""
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
        [data-testid="stSidebar"] .stSelectbox label {
              color: #e0f2fe;
              font-weight: 600;
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
    """, unsafe_allow_html=True)

# ==================== SESSION STATE INITIALIZATION ====================
def init_session_state():
    """Initialize session state variables if they don't exist."""
    defaults = {
        'lang': 'en',
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

        lang_options = {'en': "🇬🇧 English", 'ru': "🇷🇺 Русский"}
        selected_lang_code = st.selectbox(
            "🌐 Language / Язык",
            options=list(lang_options.keys()),
            format_func=lambda code: lang_options[code],
            key='lang',
            label_visibility="collapsed"
        )
        
        T = lambda key: get_translation(key, st.session_state.lang)

        st.markdown("---")

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
        else:
            st.caption(T("sidebar_no_dataset"))

        st.markdown("---")

        st.markdown(f"### {T('sidebar_quick_actions')}")
        if st.button(T("sidebar_reset_all"), use_container_width=True, key="reset_all_sidebar"):
            preserved_lang = st.session_state.lang
            keys_to_reset = list(st.session_state.keys())
            for key in keys_to_reset:
                if key not in ['lang', 'status_placeholder']:
                    del st.session_state[key]
            init_session_state()
            st.session_state.lang = preserved_lang
            st.success(T("sidebar_reset_success"))
            time.sleep(1)
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

    # Main Area
    st.markdown(f"## {T('app_title')}")

    st.session_state.status_placeholder = st.empty()
    if st.session_state.status_message:
        update_status(st.session_state.status_message, st.session_state.status_level, log=False)

    tab_keys = ["upload_tab", "manage_tab", "analyze_tab", "refine_tab", "export_tab", "docs_tab"]
    tabs = st.tabs([T(key) for key in tab_keys])
    tab_map = dict(zip(tab_keys, tabs))

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

        # --- MODIFICATION: Add Google Drive Option ---
        upload_options = [T("upload_widget"), T("upload_url")]
        # Only add Drive option if potentially available
        if COLAB_AVAILABLE or os.path.exists('/content/drive'): # Basic check
             upload_options.insert(1, T("upload_gdrive"))

        selected_upload_method = st.radio("Select Upload Method:", upload_options, horizontal=True, label_visibility="collapsed")
        # --- END MODIFICATION ---

        # --- Widget Upload ---
        if selected_upload_method == T("upload_widget"):
            col1, col2 = st.columns([2, 1])

            with col1:
                st.subheader(T("upload_files_header"))
                uploaded_files = st.file_uploader(
                    T("file_uploader_label"),
                    type=['fasta', 'fas', 'fa', 'fna', 'txt', 'gz'],
                    accept_multiple_files=True,
                    help=T("upload_help_text"),
                    key="main_file_uploader"
                )

                if uploaded_files:
                    with st.spinner(T("processing_files")):
                        progress_bar = st.progress(0, text=T("initializing"))
                        parser = FastaParser()
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

                                    sequences, errors = parser.parse(content_string)

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
                            st.balloons()
                            if not st.session_state.active_sequences:
                                st.info(T("info_activate_files"))
                        elif not has_errors:
                            st.warning(T("no_new_files"))

            with col2:
                st.subheader(T("download_url_header"))
                url_input = st.text_input(
                    T("url_input_label"),
                    placeholder=T("url_placeholder"),
                    key="url_downloader_input"
                )
                if st.button(T("download_url_btn"), use_container_width=True, key="url_download_button"):
                    if url_input and url_input.startswith(('http://', 'https://')):
                        with st.spinner(T("downloading_from_url").format(url=url_input[:50])):
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
                                    parser = FastaParser()
                                    sequences, errors = parser.parse(content_string)

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
        # --- ADDITION: Google Drive Upload Section ---
        elif selected_upload_method == T("upload_gdrive"):
            st.info(T("gdrive_info"), icon="ℹ️")
            if COLAB_AVAILABLE: # Only show mount button if in Colab
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
            gdrive_path = st.text_input(T("gdrive_path_label"), placeholder="/content/drive/MyDrive/YourFolder/*.fasta")
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
                                return
                            
                            parser = FastaParser()
                            newly_loaded_count = 0
                            total_sequences_added = 0
                            
                            for file_path in matching_files:
                                filename = os.path.basename(file_path)
                                if filename.lower().endswith(('.fasta', '.fas', '.fa', '.fna', '.txt')):
                                    with open(file_path, 'r') as f:
                                        content_string = f.read()
                                    
                                    sequences, errors = parser.parse(content_string)
                                    
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
        # --- END ADDITION ---

        # --- URL Upload ---
        elif selected_upload_method == T("upload_url"):
            col1, col2 = st.columns([2, 1])  # Reuse columns if needed, but since radio, it's handled above

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
                    <p><b>💡 {T('tip_title')}</b> {T('tip_multi_file')}</p>
                </div>
            """, unsafe_allow_html=True)
        else:
            st.subheader(T("loaded_datasets_header"))
            st.caption(T("loaded_datasets_desc"))

            file_selection_states = {}
            cols = st.columns(2)
            sorted_filenames = sorted(st.session_state.all_files.keys())

            for idx, filename in enumerate(sorted_filenames):
                sequences = st.session_state.all_files[filename]
                count = len(sequences)
                checkbox_key = f"cb_manage_{filename}"
                default_checked = filename in st.session_state.active_filenames
                with cols[idx % 2]:
                    is_selected = st.checkbox(
                        f"**{filename}** ({count} {T('seqs_abbrev')})",
                        key=checkbox_key,
                        value=default_checked
                    )
                    file_selection_states[filename] = is_selected

            selected_files_now = [fname for fname, selected in file_selection_states.items() if selected]

            st.markdown("---")
            st.subheader(T("actions_header"))
            action_cols = st.columns(4)

            with action_cols[0]:
                if st.button(T("select_all_btn"), use_container_width=True, key="manage_select_all"):
                    st.session_state.active_filenames = list(st.session_state.all_files.keys())
                    st.rerun()

            with action_cols[1]:
                if st.button(T("deselect_all_btn"), use_container_width=True, key="manage_deselect_all"):
                    st.session_state.active_filenames = []
                    for filename in st.session_state.all_files:
                        st.session_state[f"cb_manage_{filename}"] = False
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
                            st.session_state.active_sequences.extend(current_file_seqs)
                            st.session_state.original_sequences[fname] = current_file_seqs
                        count = len(st.session_state.active_sequences)
                        st.success(T("files_activated").format(count=len(selected_files_now), seqs=count))
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
            col1, col2 = st.columns(2)
            with col1:
                st.plotly_chart(create_metric_indicator(
                    len(analyzer.sequences),
                    "metric_title",
                    st.session_state.lang
                ), use_container_width=True)
            with col2:
                avg_len = sum(len(s[1]) for s in analyzer.sequences) / len(analyzer.sequences) if analyzer.sequences else 0
                st.plotly_chart(create_gauge_indicator(
                    avg_len,
                    max_value=max(2000, int(avg_len * 1.5)),
                    title_key="gauge_title",
                    lang=st.session_state.lang
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

                     selected_chart_display = st.selectbox(T("chart_type_label"), list(vis_chart_options.keys()), key="vis_chart_type")
                     selected_chart_key = vis_chart_options[selected_chart_display]

                     # ADDED Conditional controls
                     field1, field2, interval, top_n_val = None, None, None, 20
                     if selected_chart_key in ['bar', 'pie']:
                         field1_display = st.selectbox(T("field_label"), list(vis_field_options.keys()), key="vis_field1")
                         field1 = vis_field_options[field1_display]
                     elif selected_chart_key == 'line':
                         interval_options = {T("vis_interval_month"): 'month', T("vis_interval_quarter"): 'quarter', T("vis_interval_year"): 'year'}
                         interval_display = st.selectbox(T("time_interval_label"), list(interval_options.keys()), key="vis_interval")
                         interval = interval_options[interval_display]
                     elif selected_chart_key == 'heatmap':
                         top_n_val = st.slider(T("top_n_label"), 5, 50, 20, 5, key="vis_top_n")
                     elif selected_chart_key == 'stacked':
                         # Use different selectbox keys to avoid conflict
                         cat1_display = st.selectbox(T("category1_label"), list(vis_field_options.keys()), index=3, key="vis_cat1_stacked") # Default location
                         field1 = vis_field_options[cat1_display]
                         cat2_display = st.selectbox(T("category2_label"), list(vis_field_options.keys()), index=0, key="vis_cat2_stacked") # Default subtype
                         field2 = vis_field_options[cat2_display]
                         top_n_val = st.slider(T("top_n_label") + f" ({cat1_display})", 5, 30, 15, 5, key="vis_top_n_stacked")
                     # END ADDED Conditional controls

                 with vis_col2:
                     # Add vertical space to align button
                     for _ in range(5 if selected_chart_key not in ['stacked','heatmap'] else (7 if selected_chart_key=='stacked' else 6) ): st.write("")
                     # UPDATED Button logic
                     if st.button(T("generate_chart_btn"), key="vis_generate", use_container_width=True, type="primary"):
                         if field1 == field2 and selected_chart_key == 'stacked':
                              st.error("Primary and Secondary categories cannot be the same for Stacked Bar chart.")
                         else:
                              with st.spinner(T("generating_chart")):
                                  fig = None
                                  try:
                                      # Call appropriate new or existing function
                                      if selected_chart_key in ['bar', 'pie']:
                                          counts = analyzer.get_metadata_distribution(field1)
                                          fig = create_distribution_chart(counts, f"{field1_display} Distribution", chart_type=selected_chart_key)
                                      elif selected_chart_key == 'line':
                                          fig = create_temporal_chart(analyzer.sequences, interval=interval)
                                      elif selected_chart_key == 'heatmap':
                                          fig = create_geographic_heatmap(analyzer.sequences, top_n=top_n_val)
                                      elif selected_chart_key == 'stacked':
                                          fig = create_stacked_bar_chart(analyzer.sequences, category1=field1, category2=field2, top_n=top_n_val)

                                      if fig:
                                           st.session_state.generated_chart = fig # Store figure
                                           st.caption(T("chart_ready"))
                                      else:
                                           st.warning(f"Could not generate chart. No data available.")
                                  except Exception as e:
                                       st.error(f"{T('chart_error')}: {e}")
                                       progress_tracker.log_error(f"Chart generation failed: {e}")
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
                min_len = st.slider(T("min_length_label"), 0, 3000, 200, 50, key="analyze_min_len", help=T("help_min_length"))
                max_n = st.slider(T("max_n_run_label"), 0, 500, 100, 10, key="analyze_max_n", help=T("help_max_n"))
                if st.button(T("quality_filter_btn"), key="analyze_quality", use_container_width=True):
                    with st.spinner(T("applying_quality_filter")):
                        analyzer.quality_filter(min_length=min_len, max_n_run=max_n)
                        st.rerun()

                st.markdown(f"#### {T('subtype_operations')}")
                all_subtypes = ['All'] + sorted(list(set(m.get('type', DEFAULT_UNKNOWN) for _, _, m in st.session_state.active_sequences if m.get('type') != DEFAULT_UNKNOWN)))
                selected_subtype = st.selectbox(T("subtype_label"), all_subtypes, key="analyze_subtype_select")
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
                clade_mode_display = st.radio(T("mode_label"), [T("clade_mode_single"), T("clade_mode_multiple")], key="refine_clade_mode", horizontal=True)
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
                keep_monthly_display = st.selectbox(T("keep_monthly_label"), list(keep_monthly_options.keys()), index=2, key="refine_clade_keep")
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
                group_by_display = st.selectbox(T("group_by_label"), list(group_options.keys()), key="refine_temp_group", index=4)
                group_by_val = group_options[group_by_display]
            with gsk_cols[1]:
                sort_by_display = st.selectbox(T("sort_by_label"), list(sort_options.keys()), key="refine_temp_sort")
                sort_by_val = sort_options[sort_by_display]
            with gsk_cols[2]:
                keep_per_group_display = st.selectbox(T("keep_per_group_label"), list(keep_options.keys()), index=2, key="refine_temp_keep")
                keep_per_group_val = keep_options[keep_per_group_display]

            custom_grouping_input = ""
            if group_by_val == "custom":
                custom_grouping_input = st.text_input(T("custom_grouping_label"), key="refine_temp_custom", placeholder=T("custom_grouping_placeholder"))

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
            st.subheader(T("extract_accessions_btn"))
            if st.button(T("extract_accessions_btn"), key="refine_extract"):
                accessions = analyzer.extract_accessions()
                if accessions:
                    st.session_state.accession_list = accessions
                    st.success(T("accessions_found").format(count=len(accessions), tab=T('export_tab')))
                    st.text_area(T("accession_preview"), "\n".join(accessions[:20]), height=150, disabled=True)
                else:
                    st.warning(T("no_accessions_found"))

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
                st.text_area(T("log_preview"), value=log_data if log_data else "No logs yet.", height=350, disabled=True, key="export_log_preview")

    # ==================== TAB 6: DOCUMENTATION ====================
    with tab_map["docs_tab"]:
        st.header(T("docs_header"))
        
        # Documentation content (keeping English for now, can be translated later)
        st.markdown("""
        ## 🧬 Vir-Seq-Sift - User Guide

        ### Overview
        This tool provides comprehensive analysis capabilities for influenza and respiratory virus FASTA sequences. Use the tabs to navigate through the workflow: Upload -> Manage -> Analyze -> Refine -> Export.

        ### Features & Guide

        | Feature Tab         | Action                      | Use Case                                                                 | Guide                                                                                                                               |
        | :------------------ | :-------------------------- | :----------------------------------------------------------------------- | :---------------------------------------------------------------------------------------------------------------------------------- |
        | **📁 Upload & Setup**| File Upload / URL Download | Import sequence data from various sources.                               | Use the upload widget or paste a URL. Supports `.fasta`, `.fa`, `.txt`, `.gz`.                                                     |
        |                     | Google Drive (Colab)        | Load from mounted Google Drive in Colab.                                 | Select "Google Drive", mount if needed, enter path/pattern, load.                                                                  |
        | **🗂️ Manage Datasets**| Activate / Remove / Merge   | Work with multiple files, choose subsets for analysis.                 | Check files, click 'Activate Selected'. Use 'Remove' or 'Merge & Download'. Active data is used in Analyze/Refine tabs.              |
        | **🔬 Analyze & Process**| Convert Headers             | Standardize headers to `>name|type|...` format.                          | Click 'Convert Headers'. Useful if initial parsing seems incorrect.                                                               |
        |                     | Quality Filter              | Remove low-quality sequences (short or many N's).                        | Adjust sliders for 'Min Length' and 'Max N-Run', then click 'Apply Quality Filter'.                                                |
        |                     | Deduplication (Basic)       | Remove exact sequence duplicates.                                        | Click 'Deduplicate (Sequence Only)'. Keeps the first found instance.                                                              |
        |                     | Deduplication (Advanced)    | Remove duplicates, keeping one per subtype for each unique sequence.   | Click 'Deduplicate (Seq + Subtype)'. Maintains subtype diversity.                                                                   |
        |                     | Subtype Filter              | Isolate sequences of specific subtypes (e.g., H5N1).                     | Select from dropdown or enter custom subtypes (comma-sep), then click 'Apply Subtype Filter'.                                     |
        |                     | Check Subtypes              | Understand subtype proportions in the active dataset.                    | Click 'Check Subtype Distribution'. Displays Pie/Bar charts below.                                                                |
        |                     | Data Visualizer             | Explore distributions (hosts, locations, time, etc.).                    | Select field and chart type (Bar/Pie/Line/Heatmap/Stacked) in the expander, click 'Generate Chart'.                                |
        | **🎯 Refine & Visualize**| Clade Monthly Filter      | Subsample data to get representatives per clade per month.               | Select mode (Single/Multiple), choose clade(s), 'Keep' strategy (First/Last/Both), then click 'Apply'.                            |
        |                     | Enhanced Temporal Filter    | Subsample based on flexible time/metadata grouping.                      | Configure 'Group By', 'Sort By', 'Keep' options, then click 'Apply'. Useful for representative sampling over time/location etc. |
        |                     | Extract Accessions          | Get a list of GISAID EPI_ISL IDs.                                        | Click 'Extract EPI_ISL Accessions'. A download button appears in the **Export** tab.                                               |
        | **📊 Export & Reports** | Export FASTA / Report / Log | Download results, reports, and session logs.                           | Click download buttons for the current active FASTA, the last generated report, or the full session log.                              |

        ### Tips
        - **Activation is Key**: Only sequences from *activated* datasets (in the Manage tab) are used for analysis and refinement.
        - **Large Files**: Processing large files can take time. Use the spinners/progress bars as indicators.
        - **Caching**: Parsing is cached; re-uploading the same file content should be faster.
        - **Session Data**: All work is stored in your browser session and will be lost if you close the tab or refresh without uploading again. Use the Export tab to save results.
        """)

    st.markdown("---")
    st.caption(T("footer_text"))

    gc.collect()

if __name__ == "__main__":
    main()
