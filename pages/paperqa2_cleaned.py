import streamlit as st
from paperqa import Docs, Settings, agent_query
import os
import tempfile
from dotenv import load_dotenv
from datetime import datetime
import json
import requests
import urllib.parse
import time
from requests.auth import HTTPBasicAuth, HTTPDigestAuth
from urllib.parse import urljoin, urlparse
import re
from bs4 import BeautifulSoup
import xml.etree.ElementTree as ET
import asyncio
import nest_asyncio

# Streamlitで非同期を使うための設定
nest_asyncio.apply()

# Failed DOI tracking
if 'failed_dois' not in st.session_state:
    st.session_state.failed_dois = {'elsevier': [], 'other': []}

# Streamlitのキャッシュ関数
@st.cache_data
def load_api_keys():
    """APIキーの読み込みをキャッシュ"""
    env_path = "./.env"
    
    if not os.path.exists(env_path):
        st.error(f"❌ .envファイルが見つかりません: {os.path.abspath(env_path)}")
        return {"OpenAI": False, "Anthropic/Claude": False, "Google Gemini": False}
    
    success = load_dotenv(env_path)
    if not success:
        st.error(f"❌ .envファイルの読み込みに失敗しました: {os.path.abspath(env_path)}")
    else:
        st.success(f"✅ .envファイル読み込み成功: {os.path.abspath(env_path)}")
    
    return {
        "OpenAI": os.getenv("OPENAI_API_KEY") is not None,
        "Anthropic/Claude": True,  # API key is hardcoded
        "Google Gemini": os.getenv("GEMINI_API_KEY") is not None
    }

def get_paper_directory():
    """論文保存用のディレクトリパスを取得"""
    paper_dir = os.path.join(os.getcwd(), "paperqa2", "papers")
    os.makedirs(paper_dir, exist_ok=True)
    return paper_dir

def track_failed_doi(doi, reason="Download failed"):
    """失敗したDOIを記録する"""
    # Check if it's Elsevier/ScienceDirect
    is_elsevier = False
    try:
        if 'doi.org/10.1016' in doi or 'sciencedirect.com' in doi or 'elsevier.com' in doi:
            is_elsevier = True
    except:
        pass
    
    failed_entry = {
        'doi': doi,
        'reason': reason,
        'timestamp': datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    }
    
    if is_elsevier:
        if failed_entry not in st.session_state.failed_dois['elsevier']:
            st.session_state.failed_dois['elsevier'].append(failed_entry)
    else:
        if failed_entry not in st.session_state.failed_dois['other']:
            st.session_state.failed_dois['other'].append(failed_entry)

def get_comprehensive_pdf_url(paper):
    """包括的PDF URL取得（複数の方法を試行）"""
    urls_to_try = []
    
    # 1. 直接PDF URL（最優先）
    if paper.get('pdf_url'):
        urls_to_try.append(('direct_pdf', paper['pdf_url']))
    
    # 2. PMC PDF URL - 複数のパターンを試行
    if paper.get('pmcid'):
        pmcid = paper['pmcid']
        pmc_urls = [
            f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf/",
            f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf/{pmcid}.pdf"
        ]
        for pmc_url in pmc_urls:
            urls_to_try.append(('pmc_pdf', pmc_url))
    
    # 3. DOI ベース取得
    if paper.get('doi'):
        doi_url = f"https://doi.org/{paper['doi']}"
        urls_to_try.append(('doi_redirect', doi_url))
    
    return urls_to_try

def download_paper_comprehensive(paper):
    """複数のURLを試行してPDFをダウンロード"""
    urls_to_try = get_comprehensive_pdf_url(paper)
    
    for method, url in urls_to_try:
        try:
            st.write(f"🔄 Trying {method}: {url[:50]}...")
            
            # 基本的なリクエスト処理
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
            }
            
            if method == 'direct_pdf' or method == 'pmc_pdf':
                # 直接PDFダウンロード
                response = requests.get(url, headers=headers, timeout=30, allow_redirects=True)
                if response.content.startswith(b'%PDF'):
                    st.success(f"✅ {method} success!")
                    return response.content
                else:
                    st.warning(f"⚠️ {method} failed - not PDF content")
            
            elif method == 'doi_redirect':
                # DOIリダイレクトを追跡してPublisher siteへ
                pdf_content = download_from_publisher_site(url)
                if pdf_content:
                    st.success(f"✅ {method} success!")
                    return pdf_content
                else:
                    st.warning(f"⚠️ {method} failed")
                    
        except Exception as e:
            st.warning(f"⚠️ {method} failed: {str(e)[:100]}...")
            continue
    
    # All methods failed - track the DOI
    doi = paper.get('doi', url if 'doi.org' in url else 'Unknown DOI')
    track_failed_doi(doi, "All download methods failed")
    return None

def download_from_publisher_site(publisher_url):
    """Publisher siteからPDFを抽出（簡素化版）"""
    try:
        session = requests.Session()
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
            'Accept-Language': 'en-US,en;q=0.9'
        }
        
        # DOIリダイレクトの場合
        if 'doi.org' in publisher_url:
            response = session.get(publisher_url, headers=headers, allow_redirects=True, timeout=15)
            final_url = response.url
            st.write(f"🔄 Redirected to: {final_url[:50]}...")
            
            # Skip Elsevier/ScienceDirect special handling - treat as generic
            if 'sciencedirect.com' in final_url or 'elsevier.com' in final_url:
                st.write("🔄 Elsevier site detected - using generic approach...")
                track_failed_doi(publisher_url, "Elsevier site - special handling removed")
                return None
        else:
            response = session.get(publisher_url, timeout=15, headers=headers)
        
        response.raise_for_status()
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # Generic PDF link detection
        pdf_links = []
        for link in soup.find_all('a', href=True):
            href = link.get('href')
            text = link.get_text().strip().lower()
            
            if ('pdf' in href.lower() or 
                ('pdf' in text and ('download' in text or 'full' in text)) or
                ('download' in text and 'article' in text)):
                if href.startswith('/'):
                    href = urljoin(publisher_url, href)
                elif not href.startswith('http'):
                    domain = urlparse(publisher_url).netloc
                    href = urljoin(f"https://{domain}", href)
                pdf_links.append(href)
        
        # Try downloading PDFs
        for pdf_url in pdf_links[:3]:  # Try first 3 links
            try:
                st.write(f"🔄 Trying PDF: {pdf_url[:50]}...")
                pdf_response = session.get(pdf_url, timeout=30, headers=headers, allow_redirects=True)
                
                if (pdf_response.content and 
                    len(pdf_response.content) > 1000 and  
                    pdf_response.content.startswith(b'%PDF')):
                    return pdf_response.content
                    
            except Exception as e:
                st.write(f"⚠️ PDF download failed: {str(e)[:50]}...")
        
        # No valid PDF found
        track_failed_doi(publisher_url, "No valid PDF found on publisher site")
        return None
        
    except Exception as e:
        st.error(f"❌ Publisher site error: {str(e)[:100]}...")
        track_failed_doi(publisher_url, f"Publisher site error: {str(e)[:50]}")
        return None

def download_pmc_pdf(pmcid):
    """PMCからPDFをダウンロード（複数パターン試行）"""
    pmc_patterns = [
        f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf/",
        f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf/{pmcid}.pdf",
        f"https://europepmc.org/articles/{pmcid}?pdf=render"
    ]
    
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
    }
    
    for pattern in pmc_patterns:
        try:
            st.write(f"🔄 Trying PMC PDF: {pattern}")
            response = requests.get(pattern, headers=headers, timeout=30)
            
            if response.content.startswith(b'%PDF'):
                st.success("✅ PMC PDF download successful!")
                return response.content
            else:
                st.write(f"⚠️ Not PDF content: {len(response.content)} bytes")
                
        except Exception as e:
            st.write(f"⚠️ PMC pattern failed: {str(e)[:50]}...")
    
    track_failed_doi(f"PMC:{pmcid}", "PMC PDF download failed")
    return None

def search_pubmed_papers(query_terms):
    """PubMed APIを使って論文を検索"""
    try:
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        
        # クエリを構築
        search_query = " AND ".join([f'("{term}")' for term in query_terms if term.strip()])
        
        # 検索実行
        search_params = {
            'db': 'pubmed',
            'term': search_query,
            'retmax': 50,
            'retmode': 'json'
        }
        
        search_response = requests.get(f"{base_url}esearch.fcgi", params=search_params, timeout=30)
        search_data = search_response.json()
        
        if 'esearchresult' not in search_data or not search_data['esearchresult']['idlist']:
            return []
        
        pmids = search_data['esearchresult']['idlist']
        
        # 詳細情報を取得
        fetch_params = {
            'db': 'pubmed',
            'id': ','.join(pmids),
            'retmode': 'xml'
        }
        
        fetch_response = requests.get(f"{base_url}efetch.fcgi", params=fetch_params, timeout=30)
        
        # XMLパース
        root = ET.fromstring(fetch_response.content)
        
        papers = []
        for article in root.findall('.//PubmedArticle'):
            paper_info = extract_paper_info(article)
            if paper_info:
                papers.append(paper_info)
        
        return papers
        
    except Exception as e:
        st.error(f"PubMed検索エラー: {e}")
        return []

def extract_paper_info(article):
    """XMLから論文情報を抽出"""
    try:
        # タイトル
        title_elem = article.find('.//ArticleTitle')
        title = title_elem.text if title_elem is not None else "タイトル不明"
        
        # DOI
        doi = None
        for article_id in article.findall('.//ArticleId'):
            if article_id.get('IdType') == 'doi':
                doi = article_id.text
                break
        
        # PMCID
        pmcid = None
        for article_id in article.findall('.//ArticleId'):
            if article_id.get('IdType') == 'pmc':
                pmcid = article_id.text
                break
        
        # 著者
        authors = []
        for author in article.findall('.//Author'):
            last_name = author.find('LastName')
            first_name = author.find('ForeName')
            if last_name is not None and first_name is not None:
                authors.append(f"{first_name.text} {last_name.text}")
        
        # 雑誌名
        journal_elem = article.find('.//Title')
        journal = journal_elem.text if journal_elem is not None else "雑誌名不明"
        
        # 出版年
        pub_year = None
        year_elem = article.find('.//PubDate/Year')
        if year_elem is not None:
            pub_year = year_elem.text
        
        return {
            'title': title,
            'authors': ', '.join(authors[:3]),  # 最初の3人の著者
            'journal': journal,
            'year': pub_year,
            'doi': doi,
            'pmcid': pmcid
        }
        
    except Exception as e:
        st.error(f"論文情報抽出エラー: {e}")
        return None

def create_paper_entry_display(paper, index):
    """個々の論文エントリーの表示"""
    with st.container():
        col1, col2 = st.columns([4, 1])
        
        with col1:
            st.write(f"**{paper['title']}**")
            if paper.get('authors'):
                st.write(f"*著者:* {paper['authors']}")
            if paper.get('journal') and paper.get('year'):
                st.write(f"*雑誌:* {paper['journal']} ({paper['year']})")
            
            # DOI/PMCIDの表示
            identifiers = []
            if paper.get('doi'):
                identifiers.append(f"DOI: {paper['doi']}")
            if paper.get('pmcid'):
                identifiers.append(f"PMC: {paper['pmcid']}")
            if identifiers:
                st.write(f"*識別子:* {' | '.join(identifiers)}")
        
        with col2:
            download_key = f"download_{index}"
            if st.button("📄 ダウンロード", key=download_key):
                with st.spinner(f"論文をダウンロード中..."):
                    pdf_content = download_paper_comprehensive(paper)
                    
                    if pdf_content:
                        # PDFファイルを保存
                        paper_dir = get_paper_directory()
                        safe_title = re.sub(r'[^\w\s-]', '', paper['title'][:50])
                        filename = f"{safe_title}_{int(time.time())}.pdf"
                        filepath = os.path.join(paper_dir, filename)
                        
                        with open(filepath, 'wb') as f:
                            f.write(pdf_content)
                        
                        st.success(f"✅ 論文をダウンロードしました: {filename}")
                        
                        # ダウンロードボタンを提供
                        st.download_button(
                            label="💾 PDFをダウンロード",
                            data=pdf_content,
                            file_name=filename,
                            mime="application/pdf"
                        )
                    else:
                        st.error("❌ PDFのダウンロードに失敗しました")

def display_failed_dois():
    """失敗したDOIの表示"""
    st.markdown("---")
    st.header("📊 ダウンロードに失敗した論文")
    
    elsevier_count = len(st.session_state.failed_dois['elsevier'])
    other_count = len(st.session_state.failed_dois['other'])
    
    if elsevier_count == 0 and other_count == 0:
        st.success("✅ 全ての論文が正常にダウンロードされました！")
        return
    
    st.write(f"**失敗した論文数:** Elsevier: {elsevier_count}件, その他: {other_count}件")
    
    # Elsevier failures
    if elsevier_count > 0:
        st.subheader("🔬 Elsevier/ScienceDirect")
        st.info("Elsevierは特殊な処理が必要なため、手動でアクセスしてください。")
        
        for i, failed in enumerate(st.session_state.failed_dois['elsevier']):
            with st.expander(f"Elsevier論文 {i+1}: {failed['doi'][:50]}..."):
                st.write(f"**DOI/URL:** {failed['doi']}")
                st.write(f"**失敗理由:** {failed['reason']}")
                st.write(f"**時刻:** {failed['timestamp']}")
                
                # Create clickable link
                if failed['doi'].startswith('http'):
                    link_url = failed['doi']
                elif failed['doi'].startswith('10.'):
                    link_url = f"https://doi.org/{failed['doi']}"
                else:
                    link_url = failed['doi']
                
                st.markdown(f"🔗 [論文にアクセス]({link_url})")
    
    # Other failures
    if other_count > 0:
        st.subheader("📚 その他の出版社")
        
        for i, failed in enumerate(st.session_state.failed_dois['other']):
            with st.expander(f"論文 {i+1}: {failed['doi'][:50]}..."):
                st.write(f"**DOI/URL:** {failed['doi']}")
                st.write(f"**失敗理由:** {failed['reason']}")
                st.write(f"**時刻:** {failed['timestamp']}")
                
                # Create clickable link
                if failed['doi'].startswith('http'):
                    link_url = failed['doi']
                elif failed['doi'].startswith('10.'):
                    link_url = f"https://doi.org/{failed['doi']}"
                else:
                    link_url = failed['doi']
                
                st.markdown(f"🔗 [論文にアクセス]({link_url})")
    
    # Clear failed DOIs button
    if st.button("🧹 失敗リストをクリア"):
        st.session_state.failed_dois = {'elsevier': [], 'other': []}
        st.rerun()

def main():
    st.title("📚 PaperQA2 - 学術論文検索・分析システム")
    
    # API キーステータス
    api_status = load_api_keys()
    
    with st.sidebar:
        st.subheader("🔑 API Status")
        for service, available in api_status.items():
            if available:
                st.success(f"✅ {service}")
            else:
                st.error(f"❌ {service}")
        
        st.markdown("---")
        st.subheader("📁 Papers Directory")
        paper_dir = get_paper_directory()
        st.info(f"📂 {paper_dir}")
        
        # Count papers
        try:
            paper_files = [f for f in os.listdir(paper_dir) if f.endswith('.pdf')]
            st.write(f"💾 保存済み論文: {len(paper_files)}件")
        except:
            st.write("💾 保存済み論文: 0件")
    
    # メインインターfaces
    tab1, tab2, tab3 = st.tabs(["🔍 論文検索", "💬 論文分析", "📊 統計"])
    
    with tab1:
        st.header("🔍 論文検索")
        
        # 検索クエリ入力
        st.subheader("検索条件")
        query_input = st.text_area(
            "検索キーワードを入力してください（1行につき1つのキーワード）:",
            placeholder="例:\ncancer treatment\nmachine learning\nCOVID-19",
            height=100
        )
        
        if st.button("🔍 論文を検索", type="primary"):
            if not query_input.strip():
                st.warning("検索キーワードを入力してください。")
                return
            
            query_terms = [term.strip() for term in query_input.strip().split('\n') if term.strip()]
            
            with st.spinner("PubMedで論文を検索中..."):
                papers = search_pubmed_papers(query_terms)
            
            if papers:
                st.success(f"✅ {len(papers)}件の論文が見つかりました")
                
                for i, paper in enumerate(papers):
                    create_paper_entry_display(paper, i)
                    if i < len(papers) - 1:
                        st.markdown("---")
            else:
                st.warning("該当する論文が見つかりませんでした。")
    
    with tab2:
        st.header("💬 論文分析")
        
        paper_dir = get_paper_directory()
        try:
            paper_files = [f for f in os.listdir(paper_dir) if f.endswith('.pdf')]
            
            if not paper_files:
                st.info("まず「論文検索」タブで論文をダウンロードしてください。")
            else:
                st.write(f"📚 利用可能な論文: {len(paper_files)}件")
                
                # 質問入力
                question = st.text_input(
                    "論文に関する質問を入力してください:",
                    placeholder="例: この研究の主な結果は何ですか？"
                )
                
                if st.button("💬 質問する", type="primary"):
                    if not question:
                        st.warning("質問を入力してください。")
                    else:
                        with st.spinner("論文を分析中..."):
                            try:
                                settings = Settings()
                                docs = Docs(docs_path=paper_dir)
                                
                                # 非同期関数を同期的に実行
                                loop = asyncio.new_event_loop()
                                asyncio.set_event_loop(loop)
                                result = loop.run_until_complete(
                                    agent_query(question, docs=docs, settings=settings)
                                )
                                
                                st.success("✅ 分析完了")
                                st.markdown("### 📋 回答")
                                st.write(result.answer)
                                
                                if result.references:
                                    st.markdown("### 📚 参考文献")
                                    for ref in result.references:
                                        st.write(f"- {ref}")
                                
                            except Exception as e:
                                st.error(f"分析中にエラーが発生しました: {e}")
        except Exception as e:
            st.error(f"論文ディレクトリのアクセスに失敗しました: {e}")
    
    with tab3:
        st.header("📊 統計情報")
        
        paper_dir = get_paper_directory()
        try:
            paper_files = [f for f in os.listdir(paper_dir) if f.endswith('.pdf')]
            st.metric("保存済み論文数", len(paper_files))
            
            if paper_files:
                st.subheader("📁 保存済みファイル")
                for file in paper_files[:10]:  # Show first 10 files
                    st.write(f"📄 {file}")
                
                if len(paper_files) > 10:
                    st.write(f"... および {len(paper_files) - 10} 件の追加ファイル")
        
        except Exception as e:
            st.error(f"統計情報の取得に失敗しました: {e}")
    
    # Display failed DOIs at the bottom
    display_failed_dois()

if __name__ == "__main__":
    main()