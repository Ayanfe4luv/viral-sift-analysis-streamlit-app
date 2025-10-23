    # ==================== TAB 3: ANALYZE & PROCESS ====================
    with tabs[2]:
        st.markdown("### 🔬 Analyze & Process")
        
        if not st.session_state.active_sequences:
            st.markdown("""
                <div style='background: #fee2e2; padding: 30px; border-radius: 10px; border-left: 5px solid #ef4444;'>
                    <h3 style='color: #991b1b; margin-top: 0;'>❌ No Active Dataset</h3>
                    <p style='color: #7f1d1d; margin-bottom: 0;'>
                        Please activate a dataset in the <b>Manage Datasets</b> tab first.
                    </p>
                </div>
            """, unsafe_allow_html=True)
        else:
            # Dashboard metrics
            st.markdown("#### 📊 Dataset Overview")
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.metric(
                    "🧬 Total Sequences",
                    f"{len(st.session_state.active_sequences):,}",
                    help="Number of sequences in active dataset"
                )
            
            with col2:
                avg_length = sum(len(s[1]) for s in st.session_state.active_sequences) / len(st.session_state.active_sequences)
                st.metric(
                    "📏 Avg Length",
                    f"{int(avg_length):,} bp",
                    help="Average sequence length"
                )
            
            with col3:
                total_size = sum(len(s[1]) for s in st.session_state.active_sequences) / 1024 / 1024
                st.metric(
                    "💾 Total Size",
                    f"{total_size:.2f} MB",
                    help="Combined size of all sequences"
                )
            
            with col4:
                if st.session_state.original_sequences:
                    removed = len(st.session_state.original_sequences) - len(st.session_state.active_sequences)
                    st.metric(
                        "🔻 Filtered",
                        f"{removed:,}",
                        delta=f"-{removed}" if removed > 0 else "0",
                        help="Sequences removed by filters"
                    )
            
            st.markdown("---")
            
            analyzer = SequenceAnalyzer(st.session_state.active_sequences)
            
            # Processing steps with visual cards
            st.markdown("#### 🔧 Processing Steps")
            
            # Step 1: Convert Headers
            with st.expander("1️⃣ Convert Headers to Standard Format", expanded=True):
                st.markdown("""
                    <div style='background: #e0f2fe; padding: 15px; border-radius: 8px; margin-bottom: 10px;'>
                        <p style='margin: 0;'>
                            <b>What it does:</b> Standardizes headers to pipe-delimited format<br>
                            <b>Format:</b> >IsolateName|Type|Segment|Date|ID|Clade|Host|Location
                        </p>
                    </div>
                """, unsafe_allow_html=True)
                
                if st.button("🔄 Convert Headers", key="convert_btn", use_container_width=True, type="primary"):
                    with st.spinner("🔄 Converting headers..."):
                        st.session_state.active_sequences = analyzer.convert_headers()
                        st.session_state.analysis_log.append(f"✅ Converted headers at {datetime.now().strftime('%H:%M:%S')}")
                        st.success("✅ Headers converted successfully!")
                        time.sleep(0.5)
                        st.rerun()
            
            # Step 2: Quality Filter
            with st.expander("2️⃣ Apply Quality Filter", expanded=False):
                st.markdown("""
                    <div style='background: #fef3c7; padding: 15px; border-radius: 8px; margin-bottom: 10px;'>
                        <p style='margin: 0;'>
                            <b>What it does:</b> Removes sequences below quality thresholds<br>
                            <b>Filters:</b> Minimum length & Maximum N-run
                        </p>
                    </div>
                """, unsafe_allow_html=True)
                
                col1, col2 = st.columns(2)
                with col1:
                    min_length = st.slider(
                        "📏 " + T("min_length_label"),
                        0, 2000, 200, 50,
                        help="Minimum acceptable sequence length in base pairs"
                    )
                with col2:
                    max_n_run = st.slider(
                        "🧬 " + T("max_n_run_label"),
                        0, 500, 100, import streamlit as st
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

# ==================== TRANSLATIONS ====================
TRANSLATIONS = {
    "en": {
        "app_title": "🧬 FASTA Analysis Tool",
        "upload_tab": "📁 Upload & Setup",
        "manage_tab": "🗂️ Manage Datasets",
        "analyze_tab": "🔬 Analyze & Process",
        "refine_tab": "🎯 Refine & Visualize",
        "export_tab": "📊 Export & Reports",
        "docs_tab": "📖 Documentation",
        "file_uploader_label": "Upload FASTA files",
        "url_input_label": "Download from URL",
        "download_url_btn": "Download from URL",
        "convert_headers_btn": "Convert Headers",
        "quality_filter_btn": "Apply Quality Filter",
        "deduplicate_basic_btn": "Deduplicate (Sequence Only)",
        "deduplicate_advanced_btn": "Deduplicate (Sequence + Subtype)",
        "filter_subtype_btn": "Filter by Subtype",
        "check_subtypes_btn": "Check Subtype Distribution",
        "temporal_filter_btn": "Enhanced Temporal Filter",
        "clade_filter_btn": "Apply Clade Filter",
        "extract_accessions_btn": "Extract EPI_ISL Accessions",
        "export_fasta_btn": "Export FASTA",
        "export_report_btn": "Export Report",
        "no_data_msg": "No data loaded. Please upload FASTA files first.",
        "sequences_loaded": "sequences loaded",
        "processing": "Processing...",
        "complete": "Complete!",
        "min_length_label": "Min Sequence Length",
        "max_n_run_label": "Max N-Run Length",
        "subtype_label": "Select Subtype",
        "custom_subtype_placeholder": "e.g., H5N1,H3N2",
        "group_by_label": "Group By",
        "sort_by_label": "Sort By",
        "keep_label": "Keep",
        "metric_title": "Total Sequences",
        "gauge_title": "Avg Sequence Length",
        "distribution_title": "Subtype Distribution",
        "activate_btn": "Activate Selected Files",
        "merge_btn": "Merge & Download",
        "remove_btn": "Remove Selected",
        "lang_selector": "Language"
    },
    "ru": {
        "app_title": "🧬 Инструмент Анализа FASTA",
        "upload_tab": "📁 Загрузка и Настройка",
        "manage_tab": "🗂️ Управление Наборами",
        "analyze_tab": "🔬 Анализ и Обработка",
        "refine_tab": "🎯 Уточнение и Визуализация",
        "export_tab": "📊 Экспорт и Отчеты",
        "docs_tab": "📖 Документация",
        "file_uploader_label": "Загрузить файлы FASTA",
        "url_input_label": "Скачать по URL",
        "download_url_btn": "Скачать по URL",
        "convert_headers_btn": "Конвертировать Заголовки",
        "quality_filter_btn": "Применить Фильтр Качества",
        "deduplicate_basic_btn": "Дедупликация (Только Последовательность)",
        "deduplicate_advanced_btn": "Дедупликация (Последовательность + Подтип)",
        "filter_subtype_btn": "Фильтр по Подтипу",
        "check_subtypes_btn": "Проверить Распределение Подтипов",
        "temporal_filter_btn": "Улучшенный Временной Фильтр",
        "clade_filter_btn": "Применить Фильтр Клады",
        "extract_accessions_btn": "Извлечь EPI_ISL Номера",
        "export_fasta_btn": "Экспорт FASTA",
        "export_report_btn": "Экспорт Отчета",
        "no_data_msg": "Данные не загружены. Сначала загрузите файлы FASTA.",
        "sequences_loaded": "последовательностей загружено",
        "processing": "Обработка...",
        "complete": "Завершено!",
        "min_length_label": "Мин. Длина Последовательности",
        "max_n_run_label": "Макс. Длина N-Серии",
        "subtype_label": "Выбрать Подтип",
        "custom_subtype_placeholder": "например, H5N1,H3N2",
        "group_by_label": "Группировать по",
        "sort_by_label": "Сортировать по",
        "keep_label": "Оставить",
        "metric_title": "Всего Последовательностей",
        "gauge_title": "Средняя Длина Последовательности",
        "distribution_title": "Распределение Подтипов",
        "activate_btn": "Активировать Выбранные",
        "merge_btn": "Объединить и Скачать",
        "remove_btn": "Удалить Выбранные",
        "lang_selector": "Язык"
    }
}

DATE_FORMATS = ["%Y-%m-%d", "%d.%m.%Y", "%Y/%m/%d", "%Y-%m", "%Y", "%d-%b-%Y", "%b-%d-%Y", "%Y%m%d"]

# ==================== HELPER FUNCTIONS ====================
def get_translation(key, lang="en"):
    """Get translated text for a key"""
    return TRANSLATIONS.get(lang, TRANSLATIONS["en"]).get(key, key)

def parse_date(date_str):
    """Parse various date formats"""
    if isinstance(date_str, datetime):
        return date_str
    if not date_str or 'unknown' in str(date_str).lower():
        return None
    date_str = str(date_str).strip()
    for fmt in DATE_FORMATS:
        try:
            return datetime.strptime(date_str, fmt)
        except ValueError:
            continue
    return None

def update_status(message, status_type="info"):
    """Display status message"""
    if status_type == "success":
        st.success(message)
    elif status_type == "error":
        st.error(message)
    elif status_type == "warning":
        st.warning(message)
    else:
        st.info(message)

# ==================== CORE CLASSES ====================
class FastaParser:
    """Parse FASTA files and extract metadata"""
    
    def __init__(self):
        self.known_hosts = {
            'environment', 'water', 'poultry', 'wild bird', 'waterfowl', 'avian', 
            'chicken', 'turkey', 'duck', 'goose', 'human', 'swine', 'pig', 
            'equine', 'horse', 'canine', 'dog', 'cat', 'mink', 'ferret'
        }
    
    def _extract_host_and_location(self, isolate_name):
        """Extract host and location from isolate name"""
        try:
            parts = isolate_name.split('/')
            if len(parts) >= 2:
                potential_host = parts[1].lower().replace('_', ' ')
                if potential_host in self.known_hosts:
                    host = parts[1]
                    location = parts[2] if len(parts) > 2 else "Unknown"
                else:
                    host = "Unknown"
                    location = parts[1]
                return host.capitalize(), location.capitalize()
        except:
            pass
        return "Unknown", "Unknown"
    
    def _parse_header(self, header):
        """Parse FASTA header and extract metadata"""
        clean_header = header.lstrip('>').strip()
        metadata = {
            "original_header": header,
            "isolate_name": clean_header,
            "type": "Unknown",
            "segment": "Unknown",
            "collection_date": None,
            "isolate_id": "Unknown",
            "clade": "Unknown",
            "host": "Unknown",
            "location": "Unknown"
        }
        
        if '|' in clean_header:
            parts = [p.strip() for p in clean_header.split('|')]
            metadata['isolate_name'] = parts[0]
            if len(parts) > 1: metadata['type'] = parts[1]
            if len(parts) > 2: metadata['segment'] = parts[2]
            if len(parts) > 3: metadata['collection_date'] = parse_date(parts[3])
            if len(parts) > 4: metadata['isolate_id'] = parts[4]
            if len(parts) > 5: metadata['clade'] = parts[5]
        
        host, location = self._extract_host_and_location(metadata['isolate_name'])
        if metadata['host'] == "Unknown":
            metadata['host'] = host
        if metadata['location'] == "Unknown":
            metadata['location'] = location
        
        return metadata
    
    def parse(self, file_content):
        """Parse FASTA file content"""
        sequences = []
        header = None
        seq_parts = []
        
        for line in file_content.splitlines():
            if line.startswith('>'):
                if header and seq_parts:
                    sequence = "".join(seq_parts)
                    metadata = self._parse_header(header)
                    sequences.append([header, sequence, metadata])
                header = line.strip()
                seq_parts = []
            elif line.strip():
                seq_parts.append(line.strip())
        
        if header and seq_parts:
            sequence = "".join(seq_parts)
            metadata = self._parse_header(header)
            sequences.append([header, sequence, metadata])
        
        return sequences

class SequenceAnalyzer:
    """Analyze and filter FASTA sequences"""
    
    def __init__(self, sequences):
        self.sequences = list(sequences)
    
    def convert_headers(self):
        """Convert headers to standardized pipe format"""
        converted = []
        for header, seq, data in self.sequences:
            date_obj = data.get("collection_date")
            date_str = date_obj.strftime("%Y-%m-%d") if date_obj else "Unknown-Date"
            
            parts = [
                data.get('isolate_name'),
                data.get('type'),
                data.get('segment'),
                date_str,
                data.get('isolate_id'),
                data.get('clade'),
                data.get('host'),
                data.get('location')
            ]
            header_parts = [str(p) for p in parts if p and p not in ["Unknown", ""]]
            new_header = ">" + "|".join(header_parts)
            converted.append([new_header, seq, data])
        
        return converted
    
    def quality_filter(self, min_length=200, max_n_run=100):
        """Filter by sequence quality"""
        filtered = []
        removed = []
        
        for header, seq, metadata in self.sequences:
            if len(seq) < min_length:
                removed.append(header)
                continue
            
            n_runs = re.findall(r'N+', seq.upper())
            longest_n_run = max(len(run) for run in n_runs) if n_runs else 0
            
            if longest_n_run > max_n_run:
                removed.append(header)
                continue
            
            filtered.append([header, seq, metadata])
        
        return filtered, removed
    
    def deduplicate_basic(self):
        """Remove duplicate sequences"""
        seen = set()
        unique = []
        removed = []
        
        for header, seq, metadata in self.sequences:
            if seq not in seen:
                seen.add(seq)
                unique.append([header, seq, metadata])
            else:
                removed.append(header)
        
        return unique, removed
    
    def deduplicate_advanced(self):
        """Remove duplicates while preserving subtype diversity"""
        sequence_groups = defaultdict(list)
        for header, seq, metadata in self.sequences:
            sequence_groups[seq].append([header, seq, metadata])
        
        unique = []
        removed = []
        
        for seq, seq_group in sequence_groups.items():
            if len(seq_group) == 1:
                unique.extend(seq_group)
            else:
                subtype_reps = {}
                for header, seq, metadata in seq_group:
                    subtype = metadata.get('type', 'Unknown')
                    if subtype not in subtype_reps:
                        subtype_reps[subtype] = [header, seq, metadata]
                    else:
                        removed.append(header)
                unique.extend(subtype_reps.values())
        
        return unique, removed
    
    def filter_by_subtype(self, target_subtypes):
        """Filter sequences by subtype"""
        if not target_subtypes:
            return self.sequences, []
        
        target_set = {s.strip().upper() for s in target_subtypes}
        filtered = []
        removed = []
        
        for header, seq, metadata in self.sequences:
            seq_type = str(metadata.get('type', '')).strip().upper()
            if any(target in seq_type for target in target_set):
                filtered.append([header, seq, metadata])
            else:
                removed.append(header)
        
        return filtered, removed
    
    def get_subtype_distribution(self):
        """Get subtype distribution"""
        counts = Counter()
        for header, seq, metadata in self.sequences:
            subtype = metadata.get('type', 'Unknown')
            counts[subtype] += 1
        return counts
    
    def enhanced_temporal_filter(self, group_by="location_host", 
                                  sort_by="date", order="both"):
        """Enhanced temporal diversity filter"""
        df_data = []
        for header, seq, metadata in self.sequences:
            date_val = metadata.get('collection_date')
            if isinstance(date_val, str):
                date_val = parse_date(date_val)
            
            row = {
                'header': header,
                'sequence': seq,
                'date': date_val,
                'location': metadata.get('location', 'Unknown'),
                'host': metadata.get('host', 'Unknown'),
                'clade': metadata.get('clade', 'Unknown'),
                'isolate_id': metadata.get('isolate_id', 'Unknown'),
                'month': date_val.month if date_val else None
            }
            df_data.append(row)
        
        df = pd.DataFrame(df_data)
        df = df.dropna(subset=['date']).sort_values('date')
        
        # Build group key
        if group_by == "location":
            df['group_key'] = df['location']
        elif group_by == "host":
            df['group_key'] = df['host']
        elif group_by == "clade":
            df['group_key'] = df['clade']
        elif group_by == "location_host":
            df['group_key'] = df['location'] + '_' + df['host']
        else:
            df['group_key'] = df[['location', 'host', 'month', 'clade']].astype(str).agg('_'.join, axis=1)
        
        # Sort within groups
        df = df.sort_values([sort_by if sort_by != 'date' else 'date'])
        
        # Filter by order
        filtered_indices = []
        for group, gdf in df.groupby('group_key'):
            if order == "first":
                filtered_indices.append(gdf.index[0])
            elif order == "last":
                filtered_indices.append(gdf.index[-1])
            elif order == "both":
                filtered_indices.extend([gdf.index[0], gdf.index[-1]])
        
        filtered_df = df.loc[filtered_indices].drop_duplicates()
        filtered_sequences = [[row['header'], row['sequence'], {}] 
                             for _, row in filtered_df.iterrows()]
        removed = [s[0] for s in self.sequences if s[0] not in filtered_df['header'].values]
        
        return filtered_sequences, removed
    
    def extract_accessions(self):
        """Extract accession numbers"""
        accessions = []
        for header, seq, metadata in self.sequences:
            acc = metadata.get('isolate_id', '').strip()
            if acc and acc != 'Unknown':
                accessions.append(acc)
        return accessions

# ==================== VISUALIZATION FUNCTIONS ====================
def create_metric_indicator(value, title_key, lang="en"):
    """Create a metric indicator"""
    title = get_translation(title_key, lang)
    fig = go.Figure(go.Indicator(
        mode="number",
        value=value,
        title={'text': title},
        number={'font': {'size': 40}}
    ))
    fig.update_layout(height=200)
    return fig

def create_gauge_indicator(value, max_value, title_key, lang="en"):
    """Create a gauge indicator"""
    title = get_translation(title_key, lang)
    fig = go.Figure(go.Indicator(
        mode="gauge+number",
        value=value,
        title={'text': title},
        gauge={'axis': {'range': [None, max_value]},
               'bar': {'color': "darkblue"}}
    ))
    fig.update_layout(height=250)
    return fig

def create_distribution_chart(data_dict, title_key, lang="en"):
    """Create distribution pie and bar charts"""
    if not data_dict:
        return None, None
    
    title = get_translation(title_key, lang)
    df = pd.DataFrame(list(data_dict.items()), columns=['Category', 'Count'])
    
    # Pie chart
    pie_fig = px.pie(df, values='Count', names='Category', title=f"{title} (Pie)")
    
    # Bar chart
    bar_fig = px.bar(df, x='Category', y='Count', title=f"{title} (Bar)")
    
    return pie_fig, bar_fig

# ==================== CUSTOM CSS ====================
def load_custom_css():
    """Load custom CSS for better UI"""
    st.markdown("""
    <style>
    /* Main app styling */
    .main {
        background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
    }
    
    /* Card-style containers */
    .stApp [data-testid="stVerticalBlock"] > [data-testid="stVerticalBlock"] {
        background: white;
        padding: 20px;
        border-radius: 10px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        margin-bottom: 20px;
    }
    
    /* Header styling */
    h1 {
        color: #1e3a8a;
        font-weight: 700;
        text-shadow: 2px 2px 4px rgba(0,0,0,0.1);
    }
    
    h2, h3 {
        color: #2563eb;
        border-bottom: 3px solid #3b82f6;
        padding-bottom: 10px;
        margin-top: 20px;
    }
    
    /* Button styling */
    .stButton > button {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        border: none;
        border-radius: 8px;
        padding: 10px 24px;
        font-weight: 600;
        transition: all 0.3s ease;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
    }
    
    .stButton > button:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 12px rgba(0,0,0,0.15);
    }
    
    /* Success/Info boxes */
    .stSuccess, .stInfo, .stWarning, .stError {
        border-radius: 8px;
        padding: 15px;
        margin: 10px 0;
        animation: slideIn 0.3s ease;
    }
    
    @keyframes slideIn {
        from {
            opacity: 0;
            transform: translateY(-10px);
        }
        to {
            opacity: 1;
            transform: translateY(0);
        }
    }
    
    /* File uploader styling */
    .stFileUploader {
        border: 2px dashed #3b82f6;
        border-radius: 10px;
        padding: 20px;
        background: #eff6ff;
    }
    
    /* Metric cards */
    [data-testid="stMetricValue"] {
        font-size: 2.5rem;
        font-weight: 700;
        color: #1e3a8a;
    }
    
    /* Sidebar styling */
    .css-1d391kg {
        background: linear-gradient(180deg, #1e3a8a 0%, #3b82f6 100%);
    }
    
    /* Tab styling */
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
        background: white;
        border-radius: 10px;
        padding: 10px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.1);
    }
    
    .stTabs [data-baseweb="tab"] {
        border-radius: 8px;
        padding: 12px 20px;
        font-weight: 600;
    }
    
    .stTabs [aria-selected="true"] {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
    }
    
    /* Progress indicators */
    .stProgress > div > div {
        background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
    }
    
    /* Data frame styling */
    .dataframe {
        border-radius: 8px;
        overflow: hidden;
        box-shadow: 0 2px 8px rgba(0,0,0,0.1);
    }
    </style>
    """, unsafe_allow_html=True)

# ==================== SESSION STATE INITIALIZATION ====================
def init_session_state():
    """Initialize session state variables"""
    if 'lang' not in st.session_state:
        st.session_state.lang = 'en'
    if 'all_files' not in st.session_state:
        st.session_state.all_files = {}
    if 'active_sequences' not in st.session_state:
        st.session_state.active_sequences = []
    if 'original_sequences' not in st.session_state:
        st.session_state.original_sequences = []
    if 'analysis_log' not in st.session_state:
        st.session_state.analysis_log = []
    if 'processing_step' not in st.session_state:
        st.session_state.processing_step = 0

# ==================== MAIN APP ====================
def main():
    # Set page config FIRST
    st.set_page_config(
        page_title="FastaFlow - FASTA Analysis Tool",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    # Load custom CSS
    load_custom_css()
    
    init_session_state()
    
    # Sidebar with enhanced styling
    with st.sidebar:
        st.image("https://via.placeholder.com/200x80/667eea/FFFFFF?text=FastaFlow", use_container_width=True)
        st.markdown("### 🧬 FASTA Analysis")
        
        # Language selector
        lang = st.selectbox(
            "🌐 " + get_translation("lang_selector", st.session_state.lang),
            options=['en', 'ru'],
            format_func=lambda x: "🇬🇧 English" if x == 'en' else "🇷🇺 Русский",
            key='lang'
        )
        
        st.markdown("---")
        
        # Quick stats in sidebar
        if st.session_state.all_files:
            st.metric("📁 Files Loaded", len(st.session_state.all_files))
        if st.session_state.active_sequences:
            st.metric("🧬 Active Sequences", len(st.session_state.active_sequences))
            avg_len = sum(len(s[1]) for s in st.session_state.active_sequences) / len(st.session_state.active_sequences)
            st.metric("📏 Avg Length", f"{int(avg_len)} bp")
        
        st.markdown("---")
        
        # Quick actions
        st.markdown("### ⚡ Quick Actions")
        if st.button("🔄 Reset All Data", use_container_width=True):
            for key in list(st.session_state.keys()):
                del st.session_state[key]
            st.rerun()
        
        if st.session_state.active_sequences:
            if st.button("💾 Quick Export", use_container_width=True):
                output = io.StringIO()
                for header, seq, _ in st.session_state.active_sequences:
                    output.write(f"{header}\n{seq}\n")
                st.download_button(
                    "⬇️ Download FASTA",
                    data=output.getvalue(),
                    file_name=f"quick_export_{datetime.now().strftime('%Y%m%d_%H%M%S')}.fasta",
                    mime="text/plain",
                    use_container_width=True
                )
    
    T = lambda key: get_translation(key, lang)
    
    # Main title with icon
    st.markdown("""
        <div style='text-align: center; padding: 20px;'>
            <h1 style='font-size: 3rem; margin-bottom: 10px;'>
                🧬 FastaFlow
            </h1>
            <p style='font-size: 1.2rem; color: #64748b;'>
                Professional FASTA Sequence Analysis Platform
            </p>
        </div>
    """, unsafe_allow_html=True)
    
    # Tabs
    tabs = st.tabs([
        T("upload_tab"),
        T("manage_tab"),
        T("analyze_tab"),
        T("refine_tab"),
        T("export_tab"),
        T("docs_tab")
    ])
    
    # ==================== TAB 1: UPLOAD & SETUP ====================
    with tabs[0]:
        st.markdown("### 📁 Upload & Setup")
        
        # Welcome message for first-time users
        if not st.session_state.all_files:
            st.markdown("""
                <div style='background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                            color: white; padding: 30px; border-radius: 15px; text-align: center;'>
                    <h2 style='color: white; border: none;'>👋 Welcome to FastaFlow!</h2>
                    <p style='font-size: 1.1rem; margin-top: 15px;'>
                        Upload your FASTA files below to begin analyzing influenza and respiratory virus sequences
                    </p>
                </div>
            """, unsafe_allow_html=True)
            st.markdown("<br>", unsafe_allow_html=True)
        
        # File upload section with card styling
        col1, col2 = st.columns([2, 1])
        
        with col1:
            st.markdown("#### 📤 Upload Files")
            uploaded_files = st.file_uploader(
                T("file_uploader_label"),
                type=['fasta', 'fas', 'fa', 'fna', 'txt', 'gz'],
                accept_multiple_files=True,
                help="Drag and drop FASTA files here or click to browse"
            )
            
            if uploaded_files:
                with st.spinner("🔄 Processing files..."):
                    progress_bar = st.progress(0)
                    parser = FastaParser()
                    
                    for idx, uploaded_file in enumerate(uploaded_files):
                        if uploaded_file.name not in st.session_state.all_files:
                            content = uploaded_file.read()
                            if uploaded_file.name.endswith('.gz'):
                                content = gzip.decompress(content)
                            content = content.decode('utf-8')
                            sequences = parser.parse(content)
                            st.session_state.all_files[uploaded_file.name] = sequences
                        
                        progress_bar.progress((idx + 1) / len(uploaded_files))
                    
                    progress_bar.empty()
                    st.success(f"✅ Successfully loaded {len(uploaded_files)} file(s) with {sum(len(seqs) for seqs in st.session_state.all_files.values())} sequences")
                    st.balloons()
        
        with col2:
            st.markdown("#### 📊 File Statistics")
            if st.session_state.all_files:
                total_seqs = sum(len(seqs) for seqs in st.session_state.all_files.values())
                total_size = sum(sum(len(s[1]) for s in seqs) for seqs in st.session_state.all_files.values())
                
                st.metric("Files", len(st.session_state.all_files), help="Total number of loaded files")
                st.metric("Sequences", f"{total_seqs:,}", help="Total sequences across all files")
                st.metric("Total Size", f"{total_size/1024/1024:.2f} MB", help="Combined size of all sequences")
            else:
                st.info("📊 Statistics will appear here after uploading files")
        
        st.markdown("---")
        
        # URL download section
        with st.expander("🌐 Download from URL", expanded=False):
            st.markdown("**Download FASTA files directly from the web**")
            col1, col2 = st.columns([3, 1])
            
            with col1:
                url_input = st.text_input(
                    "Enter URL", 
                    placeholder="https://example.com/sequences.fasta",
                    label_visibility="collapsed"
                )
            
            with col2:
                download_btn = st.button("⬇️ Download", use_container_width=True)
            
            if download_btn and url_input:
                try:
                    with st.spinner("⏳ Downloading..."):
                        response = requests.get(url_input, timeout=30)
                        response.raise_for_status()
                        content = response.content
                        
                        filename = url_input.split('/')[-1] or 'downloaded.fasta'
                        if filename.endswith('.gz'):
                            content = gzip.decompress(content)
                        content = content.decode('utf-8')
                        
                        parser = FastaParser()
                        sequences = parser.parse(content)
                        st.session_state.all_files[filename] = sequences
                        st.success(f"✅ Downloaded {filename} ({len(sequences)} sequences)")
                        st.rerun()
                except Exception as e:
                    st.error(f"❌ Download failed: {str(e)}")
    
    # ==================== TAB 2: MANAGE DATASETS ====================
    with tabs[1]:
        st.markdown("### 🗂️ Manage Datasets")
        
        if not st.session_state.all_files:
            st.markdown("""
                <div style='background: #fef3c7; padding: 30px; border-radius: 10px; border-left: 5px solid #f59e0b;'>
                    <h3 style='color: #92400e; margin-top: 0;'>⚠️ No Files Loaded</h3>
                    <p style='color: #78350f; margin-bottom: 0;'>
                        Please upload FASTA files in the <b>Upload & Setup</b> tab first.
                    </p>
                </div>
            """, unsafe_allow_html=True)
        else:
            # File selection cards
            st.markdown("#### 📋 Select Files to Work With")
            
            selected_files = []
            
            # Create a grid layout for file cards
            cols = st.columns(2)
            for idx, (filename, sequences) in enumerate(st.session_state.all_files.items()):
                with cols[idx % 2]:
                    # Calculate file stats
                    total_length = sum(len(s[1]) for s in sequences)
                    size_mb = total_length / 1024 / 1024
                    
                    # Create card with checkbox
                    card_color = "#e0f2fe" if idx % 2 == 0 else "#f0fdf4"
                    selected = st.checkbox(
                        f"**{filename}**",
                        key=f"select_{filename}",
                        help=f"Sequences: {len(sequences):,} | Size: {size_mb:.2f} MB"
                    )
                    
                    if selected:
                        selected_files.append(filename)
                    
                    # Display file info
                    st.markdown(f"""
                        <div style='background: {card_color}; padding: 10px; border-radius: 5px; margin-bottom: 10px;'>
                            <small>
                                📊 <b>{len(sequences):,}</b> sequences | 
                                📏 <b>{size_mb:.2f}</b> MB
                            </small>
                        </div>
                    """, unsafe_allow_html=True)
            
            st.markdown("---")
            
            # Action buttons
            st.markdown("#### ⚡ Actions")
            col1, col2, col3 = st.columns(3)
            
            with col1:
                if st.button("✅ " + T("activate_btn"), use_container_width=True, type="primary"):
                    if selected_files:
                        st.session_state.active_sequences = []
                        for filename in selected_files:
                            st.session_state.active_sequences.extend(st.session_state.all_files[filename])
                        st.session_state.original_sequences = st.session_state.active_sequences.copy()
                        st.success(f"✅ Activated {len(selected_files)} file(s) with {len(st.session_state.active_sequences):,} sequences")
                        st.balloons()
                    else:
                        st.warning("⚠️ Please select at least one file")
            
            with col2:
                if st.button("🔗 " + T("merge_btn"), use_container_width=True):
                    if selected_files:
                        merged_sequences = []
                        for filename in selected_files:
                            merged_sequences.extend(st.session_state.all_files[filename])
                        
                        output = io.StringIO()
                        for header, seq, _ in merged_sequences:
                            output.write(f"{header}\n{seq}\n")
                        
                        st.download_button(
                            label=f"⬇️ Download Merged ({len(merged_sequences):,} seqs)",
                            data=output.getvalue(),
                            file_name=f"merged_{datetime.now().strftime('%Y%m%d_%H%M%S')}.fasta",
                            mime="text/plain",
                            use_container_width=True
                        )
                    else:
                        st.warning("⚠️ Please select files to merge")
            
            with col3:
                if st.button("🗑️ " + T("remove_btn"), use_container_width=True):
                    if selected_files:
                        for filename in selected_files:
                            del st.session_state.all_files[filename]
                        st.success(f"🗑️ Removed {len(selected_files)} file(s)")
                        st.rerun()
                    else:
                        st.warning("⚠️ Please select files to remove")
    
    # ==================== TAB 3: ANALYZE & PROCESS ====================
    with tabs[2]:
        st.header(T("analyze_tab"))
        
        if not st.session_state.active_sequences:
            st.warning(T("no_data_msg"))
        else:
            # Display metrics
            col1, col2 = st.columns(2)
            with col1:
                metric_fig = create_metric_indicator(
                    len(st.session_state.active_sequences),
                    "metric_title",
                    lang
                )
                st.plotly_chart(metric_fig, use_container_width=True)
            
            with col2:
                avg_length = sum(len(s[1]) for s in st.session_state.active_sequences) / len(st.session_state.active_sequences)
                gauge_fig = create_gauge_indicator(
                    avg_length,
                    5000,
                    "gauge_title",
                    lang
                )
                st.plotly_chart(gauge_fig, use_container_width=True)
            
            analyzer = SequenceAnalyzer(st.session_state.active_sequences)
            
            # Convert Headers
            if st.button(T("convert_headers_btn")):
                with st.spinner(T("processing")):
                    st.session_state.active_sequences = analyzer.convert_headers()
                    st.session_state.analysis_log.append(f"Converted headers at {datetime.now()}")
                    update_status(T("complete"), "success")
            
            # Quality Filter
            st.subheader("Quality Filter")
            col1, col2 = st.columns(2)
            with col1:
                min_length = st.slider(T("min_length_label"), 0, 2000, 200, 50)
            with col2:
                max_n_run = st.slider(T("max_n_run_label"), 0, 500, 100, 10)
            
            if st.button(T("quality_filter_btn")):
                with st.spinner(T("processing")):
                    filtered, removed = analyzer.quality_filter(min_length, max_n_run)
                    st.session_state.active_sequences = filtered
                    st.session_state.analysis_log.append(f"Quality filter: {len(removed)} removed")
                    update_status(f"Removed {len(removed)} sequences", "success")
            
            # Deduplication
            st.subheader("Deduplication")
            col1, col2 = st.columns(2)
            with col1:
                if st.button(T("deduplicate_basic_btn")):
                    with st.spinner(T("processing")):
                        unique, removed = analyzer.deduplicate_basic()
                        st.session_state.active_sequences = unique
                        st.session_state.analysis_log.append(f"Basic dedup: {len(removed)} removed")
                        update_status(f"Removed {len(removed)} duplicates", "success")
            
            with col2:
                if st.button(T("deduplicate_advanced_btn")):
                    with st.spinner(T("processing")):
                        unique, removed = analyzer.deduplicate_advanced()
                        st.session_state.active_sequences = unique
                        st.session_state.analysis_log.append(f"Advanced dedup: {len(removed)} removed")
                        update_status(f"Removed {len(removed)} duplicates", "success")
            
            # Subtype Filter
            st.subheader("Subtype Filter")
            subtype_input = st.text_input(T("subtype_label"), placeholder=T("custom_subtype_placeholder"))
            
            col1, col2 = st.columns(2)
            with col1:
                if st.button(T("filter_subtype_btn")):
                    if subtype_input:
                        subtypes = [s.strip() for s in subtype_input.split(',')]
                        with st.spinner(T("processing")):
                            filtered, removed = analyzer.filter_by_subtype(subtypes)
                            st.session_state.active_sequences = filtered
                            st.session_state.analysis_log.append(f"Subtype filter: {len(removed)} removed")
                            update_status(f"Kept {len(filtered)} sequences", "success")
            
            with col2:
                if st.button(T("check_subtypes_btn")):
                    distribution = analyzer.get_subtype_distribution()
                    pie_fig, bar_fig = create_distribution_chart(dict(distribution), "distribution_title", lang)
                    st.plotly_chart(pie_fig, use_container_width=True)
                    st.plotly_chart(bar_fig, use_container_width=True)
    
    # ==================== TAB 4: REFINE & VISUALIZE ====================
    with tabs[3]:
        st.header(T("refine_tab"))
        
        if not st.session_state.active_sequences:
            st.warning(T("no_data_msg"))
        else:
            analyzer = SequenceAnalyzer(st.session_state.active_sequences)
            
            # Temporal Filter
            st.subheader(T("temporal_filter_btn"))
            col1, col2, col3 = st.columns(3)
            
            with col1:
                group_by = st.selectbox(T("group_by_label"), 
                    ["location", "host", "clade", "location_host", "location_host_month_clade"])
            
            with col2:
                sort_by = st.selectbox(T("sort_by_label"), 
                    ["date", "location", "host", "clade"])
            
            with col3:
                order = st.selectbox(T("keep_label"), 
                    ["both", "first", "last"])
            
            if st.button(T("temporal_filter_btn")):
                with st.spinner(T("processing")):
                    filtered, removed = analyzer.enhanced_temporal_filter(group_by, sort_by, order)
                    st.session_state.active_sequences = filtered
                    st.session_state.analysis_log.append(f"Temporal filter: {len(removed)} removed")
                    update_status(f"Kept {len(filtered)} sequences", "success")
            
            # Extract Accessions
            st.subheader("Extract Accession Numbers")
            if st.button(T("extract_accessions_btn")):
                accessions = analyzer.extract_accessions()
                if accessions:
                    output = "\n".join(accessions)
                    st.download_button(
                        label=f"Download {len(accessions)} Accession Numbers",
                        data=output,
                        file_name=f"accessions_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt",
                        mime="text/plain"
                    )
                    st.success(f"Found {len(accessions)} accession numbers")
                else:
                    st.warning("No accession numbers found")
    
    # ==================== TAB 5: EXPORT & REPORTS ====================
    with tabs[4]:
        st.header(T("export_tab"))
        
        if not st.session_state.active_sequences:
            st.warning(T("no_data_msg"))
        else:
            # Export FASTA
            st.subheader("Export Processed FASTA")
            
            output = io.StringIO()
            for header, seq, _ in st.session_state.active_sequences:
                output.write(f"{header}\n{seq}\n")
            
            st.download_button(
                label=T("export_fasta_btn"),
                data=output.getvalue(),
                file_name=f"processed_{datetime.now().strftime('%Y%m%d_%H%M%S')}.fasta",
                mime="text/plain"
            )
            
            # Export Report
            st.subheader("Analysis Report")
            
            report_lines = [
                "="*70,
                "FASTA Analysis Report",
                f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
                "="*70,
                f"\nTotal Sequences: {len(st.session_state.active_sequences)}",
                f"Original Sequences: {len(st.session_state.original_sequences)}",
                f"Sequences Removed: {len(st.session_state.original_sequences) - len(st.session_state.active_sequences)}",
                "\n" + "="*70,
                "Analysis Log:",
                "="*70
            ]
            
            report_lines.extend(st.session_state.analysis_log)
            
            # Metadata summary
            analyzer = SequenceAnalyzer(st.session_state.active_sequences)
            subtype_dist = analyzer.get_subtype_distribution()
            
            report_lines.extend([
                "\n" + "="*70,
                "Subtype Distribution:",
                "="*70
            ])
            
            for subtype, count in subtype_dist.most_common():
                percentage = (count / len(st.session_state.active_sequences) * 100)
                report_lines.append(f"{subtype}: {count} ({percentage:.1f}%)")
            
            report_text = "\n".join(report_lines)
            
            st.text_area("Report Preview", report_text, height=400)
            
            st.download_button(
                label=T("export_report_btn"),
                data=report_text,
                file_name=f"analysis_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt",
                mime="text/plain"
            )
    
    # ==================== TAB 6: DOCUMENTATION ====================
    with tabs[5]:
        st.header(T("docs_tab"))
        
        st.markdown("""
        ## 🧬 FASTA Analysis Tool - User Guide
        
        ### Overview
        This tool provides comprehensive analysis capabilities for influenza and respiratory virus FASTA sequences.
        
        ### Features
        
        #### 1. Upload & Setup
        - **File Upload**: Upload multiple FASTA files (.fasta, .fas, .fa, .fna, .txt, .gz)
        - **URL Download**: Download FASTA files directly from URLs
        - Supports gzip-compressed files
        
        #### 2. Manage Datasets
        - **Multi-file Management**: Load and manage multiple FASTA files simultaneously
        - **Selective Activation**: Choose which files to analyze
        - **Merge Files**: Combine multiple files into a single dataset
        - **Remove Files**: Remove unwanted files from session
        
        #### 3. Analyze & Process
        - **Header Conversion**: Standardize headers to pipe-delimited format
        - **Quality Filtering**: Filter by minimum sequence length and maximum N-run length
        - **Deduplication**: 
          - Basic: Remove sequences with identical sequence content
          - Advanced: Preserve subtype diversity while removing duplicates
        - **Subtype Filtering**: Filter sequences by specific subtypes (e.g., H5N1, H1N1)
        - **Subtype Distribution**: Visualize subtype distribution with interactive charts
        
        #### 4. Refine & Visualize
        - **Temporal Diversity Filter**: Sample sequences based on time and metadata
          - Group by: location, host, clade, or combinations
          - Sort by: date, location, host, or clade
          - Keep: first, last, or both per group
        - **Accession Extraction**: Extract EPI_ISL and other accession numbers
        
        #### 5. Export & Reports
        - **FASTA Export**: Download processed sequences
        - **Analysis Reports**: Generate comprehensive reports with:
          - Sequence counts
          - Filtering statistics
          - Subtype distribution
          - Processing log
        
        ### Use Cases
        
        | Use Case | Steps |
        |----------|-------|
        | **Basic Quality Control** | 1. Upload files<br>2. Activate dataset<br>3. Apply quality filter<br>4. Export FASTA |
        | **Subtype Analysis** | 1. Upload files<br>2. Check subtype distribution<br>3. Filter by target subtypes<br>4. Export results |
        | **Temporal Sampling** | 1. Upload files<br>2. Configure temporal filter<br>3. Apply filter<br>4. Review and export |
        | **Deduplication** | 1. Upload files<br>2. Choose deduplication method<br>3. Apply filter<br>4. Export unique sequences |
        | **Multi-file Merge** | 1. Upload multiple files<br>2. Select files to merge<br>3. Click Merge & Download<br>4. Process merged file |
        
        ### Tips
        
        - **File Formats**: Supported formats include .fasta, .fas, .fa, .fna, .txt, and gzip-compressed versions
        - **Large Files**: The tool can handle large files, but processing may take time
        - **Metadata Parsing**: Headers are automatically parsed to extract:
          - Isolate name
          - Subtype/Type
          - Segment
          - Collection date
          - Host species
          - Location
          - Clade
          - Accession ID
        - **Language Support**: Switch between English and Russian using the sidebar dropdown
        - **Session State**: All data is stored in session state and will be lost when you close the browser
        
        ### Metadata Format
        
        The tool recognizes two main header formats:
        
        1. **Pipe-delimited**: `>IsolateName|Type|Segment|Date|ID|Clade|Host|Location`
        2. **Slash-delimited**: `>A/Host/Location/Year` (GISAID-style)
        
        ### Quality Filter Parameters
        
        - **Min Sequence Length**: Minimum acceptable sequence length (default: 200 bp)
        - **Max N-Run Length**: Maximum consecutive N bases allowed (default: 100)
        
        ### Temporal Filter Options
        
        - **Group By**: How to categorize sequences before filtering
          - `location`: Group by geographic location
          - `host`: Group by host species
          - `clade`: Group by viral clade
          - `location_host`: Group by location and host combination
          - `location_host_month_clade`: Group by all metadata fields
        
        - **Sort By**: Field used for ordering within groups
          - `date`: Sort by collection date
          - `location`, `host`, `clade`: Sort alphabetically
        
        - **Keep**: Which sequences to retain per group
          - `first`: Keep earliest/first sequence
          - `last`: Keep latest/last sequence
          - `both`: Keep both earliest and latest
        
        ### Troubleshooting
        
        - **No data loaded**: Make sure to upload files in Tab 1 and activate them in Tab 2
        - **Empty results**: Check if filters are too strict (adjust min length, max N-run)
        - **Parsing errors**: Ensure headers follow supported formats
        - **Download issues**: Use a modern browser with JavaScript enabled
        
        ### Citation
        
        If you use this tool in your research, please cite:
        
        ```
        FASTA Analysis Tool for Influenza/Respiratory Virus Research
        Version 1.0 (2025)
        ```
        
        ### Support
        
        For questions, issues, or feature requests, please contact your bioinformatics support team.
        """)

if __name__ == "__main__":
    main()
                
