import streamlit as st
import pandas as pd
import re
import duckdb
import plotly.graph_objects as go
from io import BytesIO

# ダークテーマの適用
st.set_page_config(page_title="MSP Library Search Tool", layout="wide")

# MSPファイルをパースする関数
@st.cache_data(show_spinner=False)
def parse_msp_fast(file):
    def get_match(regex, text, default=None):
        match = re.search(regex, text)
        return match.group(1) if match else default

    data = []
    content = file.read().decode("utf-8")
    records = re.split(r'\n(?=NAME: )', content)
    
    for r in records:
        if not r.strip():
            continue
        record = {}
        record['NAME'] = get_match(r"NAME: (.+)", r)
        record['PRECURSORMZ'] = float(get_match(r"PRECURSORMZ: (.+)", r))
        record['PRECURSORTYPE'] = get_match(r"PRECURSORTYPE: (.+)", r)  # PRECURSORTYPEを追加
        record['SMILES'] = get_match(r"SMILES: (.+)", r)
        record['INCHIKEY'] = get_match(r"INCHIKEY: (.+)", r)
        record['FORMULA'] = get_match(r"FORMULA: (.+)", r)
        record['RETENTIONTIME'] = float(get_match(r"RETENTIONTIME: (.+)", r))
        record['CCS'] = float(get_match(r"CCS: (.+)", r, default=0))
        record['IONMODE'] = get_match(r"IONMODE: (.+)", r)
        record['COMPOUNDCLASS'] = get_match(r"COMPOUNDCLASS: (.+)", r)
        num_peaks = int(get_match(r"Num Peaks: (\d+)", r))
        msms = re.findall(r"(\d+\.\d+)\t(\d+)", r)
        record['MSMS'] = [(float(ms), int(intensity)) for ms, intensity in msms]
        data.append(record)

    df = pd.DataFrame(data)
    return df

# 鎖長と不飽和度を抽出する関数
@st.cache_data(show_spinner=False)
def extract_chain_information(df):
    chains_info = df['NAME'].apply(lambda name: re.findall(r"(\d+):(\d+)", name))
    sn_lengths = []
    sn_unsaturations = []
    
    for chains in chains_info:
        sn1_length, sn1_unsaturation = (int(chains[0][0]), int(chains[0][1])) if len(chains) > 0 else (None, None)
        sn2_length, sn2_unsaturation = (int(chains[1][0]), int(chains[1][1])) if len(chains) > 1 else (None, None)
        sn3_length, sn3_unsaturation = (int(chains[2][0]), int(chains[2][1])) if len(chains) > 2 else (None, None)
        
        sn_lengths.append((sn1_length, sn2_length, sn3_length))
        sn_unsaturations.append((sn1_unsaturation, sn2_unsaturation, sn3_unsaturation))
    
    df['SN1_LENGTH'], df['SN2_LENGTH'], df['SN3_LENGTH'] = zip(*sn_lengths)
    df['SN1_UNSATURATION'], df['SN2_UNSATURATION'], df['SN3_UNSATURATION'] = zip(*sn_unsaturations)
    
    df['TOTAL_LENGTH'] = df[['SN1_LENGTH', 'SN2_LENGTH', 'SN3_LENGTH']].sum(axis=1, skipna=True)
    df['TOTAL_UNSATURATION'] = df[['SN1_UNSATURATION', 'SN2_UNSATURATION', 'SN3_UNSATURATION']].sum(axis=1, skipna=True)

    return df

# データクリーニング関数
@st.cache_data(show_spinner=False)
def clean_compoundclass_column(df):
    df['COMPOUNDCLASS'] = df['COMPOUNDCLASS'].str.strip().str.upper()
    return df

# データをDuckDBのテーブルにインポートし、クエリを実行する関数
@st.cache_data(show_spinner=False)
def filter_data_v3_duckdb(df, filters, mz_tolerance=0.1):
    con = duckdb.connect(database=':memory:', read_only=False)
    
    relevant_columns = ['NAME', 'COMPOUNDCLASS', 'PRECURSORMZ', 'PRECURSORTYPE', 'RETENTIONTIME',  # PRECURSORTYPEを追加
                        'SN1_LENGTH', 'SN1_UNSATURATION', 'SN2_LENGTH', 'SN2_UNSATURATION', 
                        'SN3_LENGTH', 'SN3_UNSATURATION', 'MSMS', 'SMILES']
    df_filtered = df[relevant_columns]

    con.execute("CREATE TABLE table_name AS SELECT * FROM df_filtered")

    query = """
    SELECT * FROM table_name WHERE 1=1
    """

    if filters['COMPOUNDCLASS']:
        query += f" AND COMPOUNDCLASS = '{filters['COMPOUNDCLASS'].upper()}'"
    if filters['SN1_LENGTH'] is not None:
        query += f" AND SN1_LENGTH = {filters['SN1_LENGTH']}"
    if filters['SN1_UNSATURATION'] is not None:
        query += f" AND SN1_UNSATURATION = {filters['SN1_UNSATURATION']}"
    if filters['SN2_LENGTH'] is not None:
        query += f" AND SN2_LENGTH = {filters['SN2_LENGTH']}"
    if filters['SN2_UNSATURATION'] is not None:
        query += f" AND SN2_UNSATURATION = {filters['SN2_UNSATURATION']}"
    if filters['SN3_LENGTH'] is not None:
        query += f" AND SN3_LENGTH = {filters['SN3_LENGTH']}"
    if filters['SN3_UNSATURATION'] is not None:
        query += f" AND SN3_UNSATURATION = {filters['SN3_UNSATURATION']}"
    if filters['PRECURSORMZ'] is not None:
        query += f" AND PRECURSORMZ BETWEEN {filters['PRECURSORMZ'] - mz_tolerance} AND {filters['PRECURSORMZ'] + mz_tolerance}"
    if filters['RETENTIONTIME'] is not None:
        query += f" AND RETENTIONTIME = {filters['RETENTIONTIME']}"

    filtered_df = con.execute(query).fetchdf()

    return filtered_df

# インタラクティブなMSMSスペクトルをプロットする関数（棒グラフでのスペクトル表示）
def plot_msms_spectrum_interactive(msms_data, compound_name):
    mz_values, intensity_values = zip(*msms_data)
    
    fig = go.Figure()

    fig.add_trace(go.Bar(
        x=mz_values,
        y=intensity_values,
        width=0.3,  # 幅を少し広げる
        marker=dict(color='white'),
        text=[f'{mz:.4f}' for mz in mz_values],
        textposition='outside',  # テキストをバーの外に配置
        hoverinfo='x+y',
    ))

    # 矢印付きの注釈を追加し、重ならないように調整
    annotations = []
    for i, mz in enumerate(mz_values):
        y_offset = 20 + i * 10  # 各注釈のYオフセットを少しずつずらす
        annotations.append(dict(
            x=mz,
            y=intensity_values[i],
            text=f"{mz:.4f}",
            showarrow=True,
            arrowhead=2,
            ax=0,  # 矢印の水平位置
            ay=-y_offset,  # 矢印の垂直位置（負の値で上に表示）
            font=dict(
                size=14,  # フォントサイズを大きく
                color="white"
            ),
            arrowcolor="white",
        ))

    fig.update_layout(
        title=f'MS/MS Spectrum for {compound_name}',
        xaxis_title='m/z',
        yaxis_title='Intensity',
        font=dict(
            size=16,
            color="white"
        ),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        annotations=annotations,  # 注釈を追加
        xaxis=dict(showgrid=False),
        yaxis=dict(showgrid=False),
    )

    st.plotly_chart(fig, use_container_width=True)

# SMILESから構造式を描画する関数を無効化
def draw_structure(smiles):
    st.write("SMILES structure rendering is currently disabled.")

# Streamlitアプリのレイアウト
st.title('MSP Library Search Tool')

# 初期状態のセッションステートの設定
if 'selected_name' not in st.session_state:
    st.session_state.selected_name = None
if 'filtered_data' not in st.session_state:
    st.session_state.filtered_data = pd.DataFrame()

uploaded_file = st.file_uploader("Upload MSP File", type=["msp"])

if uploaded_file is not None:
    with st.spinner("Parsing MSP file..."):
        df = parse_msp_fast(uploaded_file)

    with st.spinner("Cleaning and processing data..."):
        df = clean_compoundclass_column(df)
        df = extract_chain_information(df)

    st.subheader('Initial Data (First 5 Rows)')
    st.dataframe(df.head(5))

    filters = {
        'COMPOUNDCLASS': st.sidebar.text_input('Compound Class', key='compoundclass'),
        'SN1_LENGTH': st.sidebar.number_input('SN1 Length', min_value=0, step=1, value=None, key='sn1_length'),
        'SN1_UNSATURATION': st.sidebar.number_input('SN1 Unsaturation', min_value=0, step=1, value=None, key='sn1_unsaturation'),
        'SN2_LENGTH': st.sidebar.number_input('SN2 Length', min_value=0, step=1, value=None, key='sn2_length'),
        'SN2_UNSATURATION': st.sidebar.number_input('SN2 Unsaturation', min_value=0, step=1, value=None, key='sn2_unsaturation'),
        'SN3_LENGTH': st.sidebar.number_input('SN3 Length', min_value=0, step=1, value=None, key='sn3_length'),
        'SN3_UNSATURATION': st.sidebar.number_input('SN3 Unsaturation', min_value=0, step=1, value=None, key='sn3_unsaturation'),
        'PRECURSORMZ': st.sidebar.number_input('Precursor MZ', min_value=0.0, step=0.1, value=None, key='precursormz'),
        'RETENTIONTIME': st.sidebar.number_input('Retention Time', min_value=0.0, step=0.1, value=None, key='retentiontime')
    }
    
    mz_tolerance = st.sidebar.number_input('MZ Tolerance', min_value=0.0, step=0.01, value=0.1, key='mz_tolerance')

    # Filter Data ボタンがクリックされたときの処理
    if st.sidebar.button('Filter Data'):
        with st.spinner("Filtering data..."):
            st.session_state.filtered_data = filter_data_v3_duckdb(df, filters, mz_tolerance)

        if not st.session_state.filtered_data.empty:
            min_rt = st.session_state.filtered_data['RETENTIONTIME'].min()
            max_rt = st.session_state.filtered_data['RETENTIONTIME'].max()
            st.write(f"Retention Time min: {min_rt:.2f} ～ Max: {max_rt:.2f}")

    # 常にフィルタリング結果を表示
    st.subheader('Filtered Data')
    if not st.session_state.filtered_data.empty:
        st.dataframe(st.session_state.filtered_data)

        # 化合物選択ボックス
        selected_name = st.selectbox(
            'Select a compound to view MSMS spectrum',
            st.session_state.filtered_data['NAME'],
            index=st.session_state.filtered_data['NAME'].tolist().index(st.session_state.selected_name)
            if st.session_state.selected_name in st.session_state.filtered_data['NAME'].tolist() else 0
        )

        st.session_state.selected_name = selected_name

        # 選択された化合物のMSMSスペクトルを表示
        if st.session_state.selected_name:
            selected_row = st.session_state.filtered_data[st.session_state.filtered_data['NAME'] == st.session_state.selected_name].iloc[0]
            msms_data = selected_row['MSMS']
            smiles = selected_row['SMILES']

            plot_msms_spectrum_interactive(msms_data, st.session_state.selected_name)
            st.subheader('Molecular Structure')
            draw_structure(smiles)  # 無効化された関数が呼ばれる

    st.sidebar.header('Download Filtered Data')
    if not st.session_state.filtered_data.empty and st.sidebar.button('Download CSV'):
        csv = st.session_state.filtered_data.to_csv(index=False)
        st.sidebar.download_button(
            label="Download CSV",
            data=csv,
            file_name='filtered_data.csv',
            mime='text/csv'
        )
