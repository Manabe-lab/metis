"""
Selenium helper functions for PDF downloads
"""

import time
import random
import requests
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.common.keys import Keys
import streamlit as st
import re
from urllib.parse import urljoin, urlparse

# Check if Selenium is available
try:
    from webdriver_manager.chrome import ChromeDriverManager
    WEBDRIVER_MANAGER_AVAILABLE = True
except ImportError:
    WEBDRIVER_MANAGER_AVAILABLE = False

def human_like_click(driver, element):
    """Perform human-like click with delays and mouse movements"""
    try:
        # Scroll element into view
        driver.execute_script("arguments[0].scrollIntoView({block: 'center'});", element)
        time.sleep(random.uniform(0.5, 1.0))
        
        # Move to element with ActionChains for more realistic behavior
        actions = ActionChains(driver)
        actions.move_to_element(element)
        time.sleep(random.uniform(0.3, 0.7))  # Pause before click
        
        # Perform click
        actions.click(element)
        actions.perform()
        
        # Post-click delay
        time.sleep(random.uniform(0.5, 1.2))
        return True
    except Exception as e:
        return False

def human_like_scroll(driver):
    """Perform human-like scrolling"""
    try:
        # Random scroll amounts
        scroll_amounts = [100, 200, 300, 150, 250]
        for scroll in random.sample(scroll_amounts, 3):
            driver.execute_script(f"window.scrollTo(0, {scroll});")
            time.sleep(random.uniform(0.5, 1.0))
        
        # Scroll back to top
        driver.execute_script("window.scrollTo(0, 0);")
        time.sleep(random.uniform(0.5, 1.0))
    except:
        pass

def simulate_natural_mouse_movement(driver, start_x, start_y, end_x, end_y, duration=1.0):
    """Simulate natural mouse movement with bezier curve-like trajectory"""
    try:
        steps = max(10, int(duration * 20))  # More steps for smoother movement
        
        # Add some randomness to the path (simulate hand tremor/natural variation)
        mid_x = (start_x + end_x) / 2 + random.randint(-50, 50)
        mid_y = (start_y + end_y) / 2 + random.randint(-20, 20)
        
        for i in range(steps):
            t = i / steps
            # Quadratic bezier curve for natural movement
            x = int((1-t)**2 * start_x + 2*(1-t)*t * mid_x + t**2 * end_x)
            y = int((1-t)**2 * start_y + 2*(1-t)*t * mid_y + t**2 * end_y)
            
            driver.execute_script(f"""
                var event = new MouseEvent('mousemove', {{
                    view: window,
                    bubbles: true,
                    cancelable: true,
                    clientX: {x},
                    clientY: {y}
                }});
                document.dispatchEvent(event);
            """)
            
            # Variable speed - faster in middle, slower at ends
            speed_factor = 1 - abs(t - 0.5) * 2  # Parabolic speed curve
            sleep_time = (duration / steps) * (0.5 + speed_factor * 0.5)
            time.sleep(sleep_time + random.uniform(-0.01, 0.01))  # Small jitter
            
    except Exception as e:
        pass

def simulate_human_reading(driver, duration_range=(3, 7)):
    """Simulate human reading behavior with natural mouse movements"""
    try:
        read_time = random.uniform(*duration_range)
        st.write(f"   📖 Simulating reading behavior for {read_time:.1f}s...")
        
        # Get window dimensions for realistic movement bounds
        window_width = driver.execute_script("return window.innerWidth;")
        window_height = driver.execute_script("return window.innerHeight;")
        
        current_x = window_width // 2
        current_y = window_height // 2
        
        steps = int(read_time / 2)  # Move every ~2 seconds
        
        for step in range(steps):
            # Simulate reading - move horizontally across text, then down
            if step % 3 == 0:  # Every 3rd movement, simulate line break
                target_x = random.randint(window_width // 6, window_width // 3)  # Start of line
                target_y = current_y + random.randint(20, 40)  # Next line
            else:  # Normal reading movement
                target_x = current_x + random.randint(30, 100)  # Move right
                target_y = current_y + random.randint(-5, 5)    # Small vertical jitter
            
            # Keep within bounds
            target_x = max(50, min(target_x, window_width - 50))
            target_y = max(100, min(target_y, window_height - 100))
            
            # Natural movement to target
            movement_duration = random.uniform(0.8, 1.5)
            simulate_natural_mouse_movement(driver, current_x, current_y, target_x, target_y, movement_duration)
            
            current_x, current_y = target_x, target_y
            
            # Pause as if reading
            time.sleep(random.uniform(1.5, 2.5))
            
    except Exception as e:
        time.sleep(random.uniform(*duration_range))


def create_selenium_driver(browser_preference='auto', reuse_session=False, session_data_dir=None):
    """Create a configured Selenium driver for multiple browsers"""
    import os
    import tempfile
    import uuid
    from selenium.webdriver.firefox.options import Options as FirefoxOptions
    from selenium.webdriver.edge.options import Options as EdgeOptions
    
    # Try browsers in order: Firefox, Chrome, Edge
    browsers_to_try = ['firefox', 'chrome', 'edge'] if browser_preference == 'auto' else [browser_preference]
    
    for browser in browsers_to_try:
        try:
            st.write(f"   🌐 Attempting to use {browser.title()}...")
            
            if browser == 'firefox':
                firefox_options = FirefoxOptions()
                firefox_options.add_argument('--headless')
                firefox_options.add_argument('--width=1920')
                firefox_options.add_argument('--height=1080')
                
                # Create unique temporary profile
                temp_profile = os.path.join(tempfile.gettempdir(), f"selenium_firefox_{uuid.uuid4().hex[:8]}")
                os.makedirs(temp_profile, exist_ok=True)
                
                firefox_options.add_argument(f'--profile={temp_profile}')
                firefox_options.set_preference("dom.webdriver.enabled", False)
                firefox_options.set_preference("useAutomationExtension", False)
                
                st.write(f"   🔧 Firefox profile: {temp_profile}")
                
                if WEBDRIVER_MANAGER_AVAILABLE:
                    from webdriver_manager.firefox import GeckoDriverManager
                    from selenium.webdriver.firefox.service import Service as FirefoxService
                    service = FirefoxService(GeckoDriverManager().install())
                    driver = webdriver.Firefox(service=service, options=firefox_options)
                else:
                    driver = webdriver.Firefox(options=firefox_options)
                
                st.write(f"   ✅ Firefox driver created successfully!")
                return driver
                
            elif browser == 'chrome':
                chrome_options = Options()
                # Use non-headless mode for better compatibility with anti-bot systems
                # chrome_options.add_argument('--headless=new')  # Temporarily disable headless
                chrome_options.add_argument('--no-sandbox')
                chrome_options.add_argument('--disable-dev-shm-usage')
                chrome_options.add_argument('--disable-gpu')
                chrome_options.add_argument('--window-size=1920,1080')
                chrome_options.add_argument('--disable-blink-features=AutomationControlled')
                chrome_options.add_experimental_option("excludeSwitches", ["enable-automation"])
                chrome_options.add_experimental_option('useAutomationExtension', False)
                
                # Enhanced anti-detection measures
                chrome_options.add_argument('--user-agent=Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36')
                chrome_options.add_argument('--lang=en-US,en;q=0.9')
                chrome_options.add_argument('--accept-lang=en-US,en;q=0.9')
                chrome_options.add_argument('--disable-web-security')
                chrome_options.add_argument('--disable-extensions')
                chrome_options.add_argument('--disable-plugins')
                chrome_options.add_argument('--disable-images')  # Speed up loading
                
                # Additional stealth options
                prefs = {
                    "profile.default_content_setting_values": {
                        "notifications": 2,  # Block notifications
                        "plugins": 2,       # Block plugins
                        "popups": 2,        # Block popups
                        "geolocation": 2,   # Block location sharing
                        "media_stream": 2   # Block camera/mic
                    },
                    "profile.managed_default_content_settings": {
                        "images": 2  # Block images for faster loading
                    }
                }
                chrome_options.add_experimental_option("prefs", prefs)
                
                # Handle session reuse
                if reuse_session and session_data_dir and os.path.exists(session_data_dir):
                    # Reuse existing session data directory
                    temp_profile = session_data_dir
                    st.write(f"   🔄 Reusing session data: {temp_profile}")
                else:
                    # Create unique temporary profile
                    temp_profile = os.path.join(tempfile.gettempdir(), f"selenium_chrome_{uuid.uuid4().hex[:8]}")
                    os.makedirs(temp_profile, exist_ok=True)
                    
                    # Store the profile path for potential reuse
                    if st.session_state is not None:
                        st.session_state.chrome_session_dir = temp_profile
                
                chrome_options.add_argument(f'--user-data-dir={temp_profile}')
                chrome_options.add_argument('--profile-directory=Default')
                chrome_options.add_argument('--disable-web-security')
                chrome_options.add_argument('--disable-features=VizDisplayCompositor')
                chrome_options.add_argument('--no-first-run')
                chrome_options.add_argument('--disable-default-apps')
                
                st.write(f"   🔧 Chrome profile: {temp_profile}")
                
                if WEBDRIVER_MANAGER_AVAILABLE:
                    from webdriver_manager.chrome import ChromeDriverManager
                    from selenium.webdriver.chrome.service import Service as ChromeService
                    service = ChromeService(ChromeDriverManager().install())
                    driver = webdriver.Chrome(service=service, options=chrome_options)
                else:
                    driver = webdriver.Chrome(options=chrome_options)
                
                # Override navigator.webdriver and add stealth measures
                driver.execute_script("""
                    Object.defineProperty(navigator, 'webdriver', {get: () => undefined});
                    
                    // Override the `plugins` property to use a custom getter
                    Object.defineProperty(navigator, 'plugins', {
                        get: () => [1, 2, 3, 4, 5]  // Fake plugins
                    });
                    
                    // Override the `languages` property to use a custom getter
                    Object.defineProperty(navigator, 'languages', {
                        get: () => ['en-US', 'en']
                    });
                    
                    // Override chrome runtime
                    window.chrome = {
                        runtime: {}
                    };
                    
                    // Mock permissions
                    Object.defineProperty(navigator, 'permissions', {
                        get: () => ({
                            query: () => Promise.resolve({ state: 'granted' })
                        })
                    });
                """)
                st.write(f"   ✅ Chrome driver created successfully with stealth measures!")
                return driver
                
            elif browser == 'edge':
                edge_options = EdgeOptions()
                edge_options.add_argument('--headless')
                edge_options.add_argument('--no-sandbox')
                edge_options.add_argument('--disable-dev-shm-usage')
                edge_options.add_argument('--disable-gpu')
                edge_options.add_argument('--window-size=1920,1080')
                edge_options.add_argument('--disable-blink-features=AutomationControlled')
                edge_options.add_experimental_option("excludeSwitches", ["enable-automation"])
                edge_options.add_experimental_option('useAutomationExtension', False)
                
                # Create unique temporary profile
                temp_profile = os.path.join(tempfile.gettempdir(), f"selenium_edge_{uuid.uuid4().hex[:8]}")
                os.makedirs(temp_profile, exist_ok=True)
                
                edge_options.add_argument(f'--user-data-dir={temp_profile}')
                edge_options.add_argument('--profile-directory=Default')
                
                st.write(f"   🔧 Edge profile: {temp_profile}")
                
                if WEBDRIVER_MANAGER_AVAILABLE:
                    from webdriver_manager.microsoft import EdgeChromiumDriverManager
                    from selenium.webdriver.edge.service import Service as EdgeService
                    service = EdgeService(EdgeChromiumDriverManager().install())
                    driver = webdriver.Edge(service=service, options=edge_options)
                else:
                    driver = webdriver.Edge(options=edge_options)
                
                # Override navigator.webdriver
                driver.execute_script("Object.defineProperty(navigator, 'webdriver', {get: () => undefined})")
                st.write(f"   ✅ Edge driver created successfully!")
                return driver
                
        except Exception as browser_error:
            st.write(f"   ⚠️ {browser.title()} failed: {str(browser_error)[:50]}...")
            continue
    
    # If all browsers fail
    raise Exception("No compatible browser found. Please install Chrome, Firefox, or Edge with their respective drivers.")


def create_persistent_chrome_driver():
    """Create a Chrome driver with persistent session data for manual authentication"""
    import os
    import tempfile
    import uuid
    
    try:
        st.write("🔧 Creating persistent Chrome session for manual authentication...")
        
        # Create a persistent session directory
        session_id = str(uuid.uuid4().hex[:8])
        session_dir = os.path.join(tempfile.gettempdir(), f"manual_auth_chrome_{session_id}")
        os.makedirs(session_dir, exist_ok=True)
        
        chrome_options = Options()
        # Don't use headless for manual authentication
        chrome_options.add_argument('--no-sandbox')
        chrome_options.add_argument('--disable-dev-shm-usage')
        chrome_options.add_argument('--disable-gpu')
        chrome_options.add_argument('--window-size=1920,1080')
        
        # Use persistent profile
        chrome_options.add_argument(f'--user-data-dir={session_dir}')
        chrome_options.add_argument('--profile-directory=Default')
        
        # Create driver
        if WEBDRIVER_MANAGER_AVAILABLE:
            from webdriver_manager.chrome import ChromeDriverManager
            from selenium.webdriver.chrome.service import Service as ChromeService
            service = ChromeService(ChromeDriverManager().install())
            driver = webdriver.Chrome(service=service, options=chrome_options)
        else:
            driver = webdriver.Chrome(options=chrome_options)
        
        # Store session info
        st.session_state.manual_auth_session_dir = session_dir
        st.session_state.manual_auth_driver = driver
        
        st.success(f"✅ Chrome driver created with persistent session: {session_dir}")
        return driver, session_dir
        
    except Exception as e:
        st.error(f"❌ Failed to create persistent Chrome driver: {str(e)}")
        return None, None


def download_pmc_with_selenium(pmc_url):
    """Download PDF from PMC using Selenium to handle POW challenge"""
    driver = None
    try:
        # PMC URL処理 - より安全な実装
        original_url = pmc_url
        
        # PMCIDを抽出
        pmcid_match = re.search(r'PMC\d+', pmc_url)
        if pmcid_match:
            pmcid = pmcid_match.group()
            # 標準的なPMC URLを構築
            pmc_url = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/"
            pdf_url = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/pdf/"
        else:
            # PMCIDが見つからない場合の処理
            if '/pdf/' in pmc_url:
                pdf_url = pmc_url
                pmc_url = pmc_url.replace('/pdf/', '/')
            else:
                pmc_url = original_url
                pdf_url = f"{pmc_url.rstrip('/')}/pdf/"
        
        st.write("🤖 Using Selenium for PMC PDF download...")
        
        driver = create_selenium_driver()
        
        # First visit article page
        st.write(f"   1️⃣ Visiting article page...")
        driver.get(pmc_url)
        time.sleep(2)
        
        # Now navigate to PDF
        st.write(f"   2️⃣ Navigating to PDF URL...")
        driver.get(pdf_url)
        
        # Wait for POW challenge
        st.write(f"   3️⃣ Waiting for POW challenge to complete...")
        max_wait = 30
        wait_interval = 2
        total_waited = 0
        
        while total_waited < max_wait:
            time.sleep(wait_interval)
            total_waited += wait_interval
            
            page_source = driver.page_source
            
            # Check if POW completed
            if "Preparing to download" not in page_source:
                st.write(f"   ✅ POW challenge completed after {total_waited}s")
                break
            
            if total_waited % 6 == 0:  # Update every 6 seconds
                st.write(f"   ⏳ Still waiting... ({total_waited}/{max_wait}s)")
        
        # Get cookies
        cookies = driver.get_cookies()
        pow_cookies = [c for c in cookies if 'pow' in c['name'].lower()]
        
        if pow_cookies:
            st.write(f"   4️⃣ Found POW cookies, downloading PDF...")
            
            # Create session with cookies
            session = requests.Session()
            for cookie in cookies:
                session.cookies.set(
                    cookie['name'], 
                    cookie['value'],
                    domain=cookie.get('domain', '.ncbi.nlm.nih.gov'),
                    path=cookie.get('path', '/')
                )
            
            # Download with cookies
            user_agent = driver.execute_script("return navigator.userAgent;")
            headers = {
                'User-Agent': user_agent,
                'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
                'Accept-Language': 'en-US,en;q=0.9',
                'Referer': pmc_url
            }
            
            response = session.get(pdf_url, headers=headers, timeout=30)
            
            if response.content.startswith(b'%PDF'):
                st.success("   ✅ Successfully downloaded PDF with Selenium!")
                return response.content
            else:
                st.warning("   ⚠️ Downloaded content is not a PDF")
        else:
            st.warning("   ⚠️ No POW cookies found")
            
    except Exception as e:
        st.error(f"   ❌ Selenium error: {str(e)}")
    finally:
        if driver:
            driver.quit()
    
    return None


def download_sciencedirect_with_selenium(url):
    """Download PDF from ScienceDirect using Selenium"""
    driver = None
    try:
        st.write("🤖 Using Selenium for ScienceDirect PDF download...")
        
        driver = create_selenium_driver()
        
        # First, establish session by visiting main site with human-like behavior
        st.write(f"   0️⃣ Establishing session with ScienceDirect...")
        driver.get("https://www.sciencedirect.com")
        
        # Minimal waiting to avoid triggering Cloudflare
        st.write("   ⏳ Brief wait for initial page load...")
        time.sleep(random.uniform(2, 3))  # Reduced from human-like behavior
        
        # Accept cookies with human-like behavior
        try:
            # Wait for cookie banner to load
            wait = WebDriverWait(driver, 15)
            cookie_button = wait.until(EC.element_to_be_clickable((By.ID, "onetrust-accept-btn-handler")))
            
            # Use human-like click
            if human_like_click(driver, cookie_button):
                st.write("   🍪 Accepted cookies with human-like behavior")
            else:
                # Fallback to regular click
                cookie_button.click()
                st.write("   🍪 Accepted cookies (fallback)")
                time.sleep(random.uniform(1, 2))
        except:
            # Try alternative cookie selectors
            try:
                cookie_alternatives = [
                    "#onetrust-accept-btn-handler",
                    ".onetrust-accept-btn-handler", 
                    "[data-id='accept-all-cookies']",
                    "button[id*='cookie'][id*='accept']",
                    "button[class*='cookie'][class*='accept']"
                ]
                for selector in cookie_alternatives:
                    try:
                        cookie_btn = driver.find_element(By.CSS_SELECTOR, selector)
                        if cookie_btn.is_displayed():
                            # Use human-like click for alternatives too
                            if human_like_click(driver, cookie_btn):
                                st.write("   🍪 Accepted cookies (alternative method with human-like behavior)")
                            else:
                                cookie_btn.click()
                                st.write("   🍪 Accepted cookies (alternative method)")
                                time.sleep(random.uniform(1, 2))
                            break
                    except:
                        continue
            except:
                st.write("   ⚠️ No cookie banner found")
        
        # Additional human-like behavior after cookie handling
        time.sleep(random.uniform(1, 2))
        
        # Visit the URL (DOI or direct) with human-like delay
        st.write(f"   1️⃣ Navigating to article...")
        
        # Human-like delay before navigation
        time.sleep(random.uniform(2, 4))
        driver.get(url)
        
        # Brief pause instead of human-like reading to avoid Cloudflare
        time.sleep(random.uniform(1, 2))
        
        # Wait for page to load and handle Cloudflare with human-like behavior
        st.write(f"   ⏳ Waiting for page to fully load...")
        time.sleep(3)  # Initial wait
        
        # Check for Cloudflare challenge - offer immediate manual override
        max_cloudflare_wait = 6  # Very short wait - 6 seconds only
        cloudflare_wait = 0
        cloudflare_detected = False
        
        while cloudflare_wait < max_cloudflare_wait:
            page_source = driver.page_source.lower()
            page_title = driver.title.lower()
            current_url = driver.current_url.lower()
            
            # Check multiple security/login challenge indicators
            challenge_indicators = [
                'cloudflare' in page_source,
                'checking your browser' in page_source,
                'just a moment' in page_title,
                'please wait' in page_title,
                'security check' in page_source,
                'ddos protection' in page_source,
                'ray id' in page_source,
                # Institutional login indicators
                'institution/login' in current_url,
                'institutional access' in page_source,
                'sign in' in page_title.lower(),
                'login' in page_title.lower() and 'science' in page_source.lower(),
                'authentication required' in page_source,
                'access denied' in page_source,
                'subscription required' in page_source
            ]
            
            if any(challenge_indicators):
                cloudflare_detected = True
                # Determine the specific challenge type
                if 'institution/login' in current_url:
                    st.error(f"   🏛️ **Institutional login page detected!**")
                    challenge_type = "institutional login"
                elif any(['cloudflare' in page_source, 'just a moment' in page_title]):
                    st.error(f"   🚨 **Cloudflare challenge detected!**")
                    challenge_type = "Cloudflare challenge"
                else:
                    st.error(f"   🔒 **Access challenge detected!**")
                    challenge_type = "access challenge"
                    
                st.warning(f"   🚫 **Switching to manual mode for {challenge_type}**")
                break  # Exit the loop immediately
                
            time.sleep(1)
            cloudflare_wait += 1
        
        # Handle the result
        if cloudflare_detected:
            # Get current URL and window handle info
            current_url = driver.current_url
            window_title = driver.title
            
            # Determine challenge type for display
            if 'institution/login' in current_url.lower():
                challenge_display = "🏛️ **INSTITUTIONAL LOGIN REQUIRED**"
                challenge_desc = "institutional login"
                action_needed = "complete the institutional authentication"
            else:
                challenge_display = "🚨 **CLOUDFLARE CHALLENGE DETECTED**"  
                challenge_desc = "Cloudflare challenge"
                action_needed = "complete the Cloudflare challenge"
            
            # Try to bring browser window to front (may not work in all environments)
            try:
                driver.switch_to.window(driver.current_window_handle)
                # Try JavaScript to request focus
                driver.execute_script("window.focus();")
                # Try to maximize window to make it more visible
                driver.maximize_window()
            except Exception as e:
                st.write(f"   ⚠️ Could not automatically focus browser window: {str(e)[:50]}...")
            
            # Create a visually prominent alert
            alert_color = "#ff4444" if "CLOUDFLARE" in challenge_display else "#ff6600"
            st.markdown(f"""
            <div style="background-color: {alert_color}; color: white; padding: 20px; border-radius: 10px; margin: 10px 0; text-align: center;">
                <h2>{challenge_display}</h2>
                <h3>⚡ IMMEDIATE ACTION REQUIRED ⚡</h3>
            </div>
            """, unsafe_allow_html=True)
            
            st.markdown(f"""
            <div style="background-color: #ffeeee; padding: 15px; border-radius: 5px; border-left: 5px solid {alert_color};">
            
            **🔥 A browser window is open and needs your attention!**
            
            **📋 Browser Window Info:**
            - **Title:** `{window_title}`
            - **URL:** `{current_url[:60]}...`
            
            **👆 HOW TO SWITCH TO THE BROWSER:**
            
            </div>
            """, unsafe_allow_html=True)
            
            # Different instructions based on challenge type
            if 'institution/login' in current_url.lower():
                instructions = f"""
                **🏛️ Institutional Login Instructions:**
                
                1. **First:** Click the blue button below to prepare window switching
                2. **Then:** Use **Alt+Tab** (Windows/Linux) or **Cmd+Tab** (Mac) to find the browser window
                3. **Look for:** Browser window with institutional login page
                4. **Complete:** Enter your institutional credentials or select your institution
                5. **Wait:** Until you see the actual article content (not the login page)
                6. **Return:** To this page and click the green Continue button
                
                ⚠️ **IMPORTANT:** Complete your institutional authentication first!
                """
            else:
                instructions = """
                **🔧 Cloudflare Challenge Instructions:**
                
                1. **First:** Click the blue button below to prepare window switching
                2. **Then:** Use **Alt+Tab** (Windows/Linux) or **Cmd+Tab** (Mac) to find the browser window
                3. **Look for:** Browser window with title containing "Just a moment" or "Cloudflare"
                4. **Complete:** Click the "I'm not a robot" checkbox 
                5. **Wait:** Until you see the actual article content (not the challenge)
                6. **Return:** To this page and click the green Continue button
                
                ⚠️ **IMPORTANT:** Do NOT close the browser - keep it open!
                """
            
            st.info(instructions)
            
            # STOP ALL AUTOMATION - User must take action
            st.markdown("---")
            st.markdown("""
            <div style="background-color: #ff6666; color: white; padding: 15px; border-radius: 5px; text-align: center; font-size: 18px;">
                <strong>⏹️ AUTOMATION COMPLETELY STOPPED ⏹️</strong><br>
                Manual action required to continue
            </div>
            """, unsafe_allow_html=True)
            
            st.markdown("---")
            
            col1, col2, col3 = st.columns([1,2,1])
            with col2:
                if st.button("🔍 STEP 1: Click here, then switch to browser with Alt+Tab", 
                           key="focus_help", type="primary", use_container_width=True):
                    st.success("✅ Step 1 complete! Now use **Alt+Tab** to switch to the browser!")
                    st.session_state.step1_done = True
                    
            # Only show step 2 after step 1 is completed
            if st.session_state.get('step1_done', False):
                st.markdown("---")
                col1, col2, col3 = st.columns([1,2,1])
                with col2:
                    resume_text = f"✅ STEP 2: I completed the {challenge_desc} - Resume"
                    if st.button(resume_text, 
                               key="cloudflare_override", type="secondary", use_container_width=True):
                        st.success("🎉 Manual intervention completed! Resuming automation...")
                        # Clear the step flag
                        st.session_state.step1_done = False
                        # Brief pause to ensure page is ready
                        time.sleep(2)
                        # Continue with the function (don't return None)
                    else:
                        st.warning(f"⏸️ **Complete the {challenge_desc} in the browser, then click Step 2**")
                        # Force function to exit and wait for user
                        raise Exception("CLOUDFLARE_WAIT_USER_ACTION")
            else:
                st.error(f"❌ **Complete Step 1 first, then the {challenge_desc} in the browser window**")
                # Force function to exit and wait for user  
                raise Exception("CLOUDFLARE_WAIT_USER_ACTION")
        else:
            st.write("   ✅ Page loaded successfully, no security challenges detected")
        
        # Get current URL after redirects
        current_url = driver.current_url
        st.write(f"   📍 Current URL: {current_url[:80]}...")
        
        # Look for PDF download button/link
        st.write(f"   2️⃣ Looking for PDF download option...")
        
        pdf_downloaded = False
        pdf_url = None
        
        # Try multiple selectors for PDF links - Enhanced for ScienceDirect
        pdf_selectors = [
            "a[aria-label*='Download PDF']",
            "a[title*='Download PDF']",
            "a[title*='PDF']",
            "a.download-pdf-link",
            "button.pdf-download",
            "a[href*='pdfft']",
            "a[href*='/pdf/']",
            "button[title*='PDF']",
            "button[aria-label*='PDF']",
            "span:contains('Download PDF')",
            "span:contains('View PDF')",
            "span:contains('PDF')",
            "button:contains('PDF')",
            "button:contains('Download PDF')",
            "a:contains('Download PDF')",
            "a:contains('View PDF')",
            ".pdf-download-btn",
            "a.pdf-download",
            ".download-link",
            "[data-testid*='pdf']",
            "[data-aa-name*='pdf']",
            "a.anchor-text[href*='pdf']",
            ".icon-pdf + a",
            "button.btn-primary[title*='PDF']",
            # Additional ScienceDirect patterns for Download PDF links
            "a[href*='download'][href*='pdf']",
            "button[data-aa-name*='download']",
            "a[data-aa-name*='download-pdf']",
            "*[role='button'][aria-label*='Download']",
            "div[data-testid='download-pdf'] a",
            "div[data-testid='download-pdf'] button"
        ]
        
        for selector in pdf_selectors:
            try:
                # Wait for element
                wait = WebDriverWait(driver, 5)
                if ':contains(' in selector:
                    # Use XPath for text content
                    xpath = f"//*[contains(text(), '{selector.split('(')[1].split(')')[0]}')]"
                    pdf_element = wait.until(EC.presence_of_element_located((By.XPATH, xpath)))
                else:
                    pdf_element = wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, selector)))
                
                st.write(f"   ✅ Found PDF element with selector: {selector}")
                
                # Get href if it's a link
                if pdf_element.tag_name == 'a':
                    pdf_url = pdf_element.get_attribute('href')
                    if pdf_url:
                        st.write(f"   🔗 PDF URL: {pdf_url[:80]}...")
                        break
                else:
                    # Click button
                    st.write(f"   🖱️ Clicking PDF download button...")
                    pdf_element.click()
                    time.sleep(3)
                    pdf_downloaded = True
                    break
                    
            except:
                continue
        
        # Enhanced PDF detection for ScienceDirect
        if not pdf_url and not pdf_downloaded:
            st.write(f"   3️⃣ Enhanced ScienceDirect PDF detection...")
            
            # Method 1: Look for any download-related elements
            try:
                download_elements = driver.find_elements(By.XPATH, "//a[contains(@href, 'pdf') or contains(text(), 'PDF') or contains(text(), 'Download')]")
                if download_elements:
                    for elem in download_elements:
                        href = elem.get_attribute('href')
                        text = elem.text.lower()
                        if href and ('pdf' in href or 'download' in text):
                            pdf_url = href
                            st.write(f"   ✅ Found download element: {text} -> {href[:60]}...")
                            break
            except:
                pass
            
            # Method 2: Check for ScienceDirect specific patterns
            if not pdf_url:
                # Look for "View PDF" or similar buttons
                try:
                    view_pdf_buttons = driver.find_elements(By.XPATH, "//button[contains(text(), 'PDF') or contains(@title, 'PDF')]")
                    if view_pdf_buttons:
                        st.write(f"   🖱️ Found PDF button, attempting click...")
                        view_pdf_buttons[0].click()
                        time.sleep(2)
                        
                        # Check if new tab opened
                        if len(driver.window_handles) > 1:
                            driver.switch_to.window(driver.window_handles[-1])
                            pdf_url = driver.current_url
                            st.write(f"   ✅ PDF opened in new tab: {pdf_url[:60]}...")
                            driver.switch_to.window(driver.window_handles[0])  # Switch back
                except:
                    pass
            
            # Method 3: Construct PDF URL from article URL
            if not pdf_url:
                st.write(f"   🔧 Constructing PDF URL from article URL...")
                
                # Extract PII from current URL
                pii_match = re.search(r'/pii/([A-Z0-9]+)', current_url)
                if pii_match:
                    pii = pii_match.group(1)
                    pdf_url = f"https://www.sciencedirect.com/science/article/pii/{pii}/pdfft"
                    st.write(f"   🔗 Constructed PDF URL: {pdf_url}")
        
        # Get cookies and headers
        cookies = driver.get_cookies()
        user_agent = driver.execute_script("return navigator.userAgent;")
        
        # Create session with cookies
        session = requests.Session()
        for cookie in cookies:
            session.cookies.set(
                cookie['name'], 
                cookie['value'],
                domain=cookie.get('domain', '.sciencedirect.com'),
                path=cookie.get('path', '/')
            )
        
        # Try to download PDF
        if pdf_url:
            st.write(f"   4️⃣ Downloading PDF...")
            
            headers = {
                'User-Agent': user_agent,
                'Accept': 'application/pdf,text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8',
                'Accept-Language': 'en-US,en;q=0.9',
                'Accept-Encoding': 'gzip, deflate, br',
                'DNT': '1',
                'Connection': 'keep-alive',
                'Upgrade-Insecure-Requests': '1',
                'Sec-Fetch-Dest': 'document',
                'Sec-Fetch-Mode': 'navigate',
                'Sec-Fetch-Site': 'same-origin',
                'Cache-Control': 'max-age=0',
                'Referer': current_url
            }
            
            response = session.get(pdf_url, headers=headers, timeout=30, allow_redirects=True)
            
            if response.content.startswith(b'%PDF'):
                st.success("   ✅ Successfully downloaded PDF with Selenium!")
                return response.content
            else:
                st.warning(f"   ⚠️ Response is not a PDF (Status: {response.status_code})")
                
                # Check if it's access denied
                if response.status_code == 403 or 'access denied' in response.text.lower():
                    st.error("   🔒 Access denied - subscription required")
        else:
            st.warning("   ⚠️ Could not find PDF download link")
    
    except Exception as e:
        if "CLOUDFLARE_WAIT_USER_ACTION" in str(e):
            # This is expected - user needs to complete Cloudflare challenge
            # Don't close the driver, keep it open for user interaction
            st.info("🔄 **Browser window kept open for manual Cloudflare completion**")
            # Store driver in session state so it persists
            st.session_state.selenium_driver = driver
            return None
        else:
            st.error(f"   ❌ Selenium error: {str(e)}")
            if driver:
                driver.quit()
            return None
    finally:
        # Don't close driver if waiting for Cloudflare completion
        if 'e' in locals() and "CLOUDFLARE_WAIT_USER_ACTION" in str(e):
            st.info("🔄 **Browser session maintained for future use**")
            pass  # Keep driver open
        elif driver and not st.session_state.get('step1_done', False):
            # Only quit if we're not in a Cloudflare flow
            try:
                # Save session info before closing (optional)
                st.info("💾 **Session cookies should be preserved for next access**")
                driver.quit()
            except:
                pass
    
    return None


def download_doi_with_selenium_with_driver(url, existing_driver):
    """Download PDF using an existing pre-authenticated driver"""
    try:
        st.write("🔐 Using pre-authenticated browser session...")
        
        # Use the existing driver that should already be authenticated
        driver = existing_driver
        
        # Check current URL - if already on the target, we're good
        current_url = driver.current_url
        if url.lower() in current_url.lower() or current_url.lower() in url.lower():
            st.success("✅ Already on target page from pre-authentication!")
        else:
            st.write("📍 Navigating to target URL...")
            driver.get(url)
            time.sleep(2)
        
        # Skip the Cloudflare detection since user should have already cleared it
        st.write("⚡ Skipping authentication checks - assuming pre-cleared...")
        
        # Continue with PDF detection and download logic
        # Use the same logic as the main function but skip authentication
        current_url = driver.current_url
        st.write(f"📍 Current URL: {current_url[:80]}...")
        
        # Try ScienceDirect specific approach if it's a ScienceDirect URL
        if 'sciencedirect.com' in current_url.lower():
            return download_sciencedirect_with_selenium(current_url)
        else:
            # For other sites, implement basic PDF detection
            st.write("🔍 Searching for PDF download links...")
            
            # Look for PDF links
            try:
                pdf_links = driver.find_elements(By.PARTIAL_LINK_TEXT, "PDF")
                if not pdf_links:
                    pdf_links = driver.find_elements(By.CSS_SELECTOR, "a[href*='pdf']")
                
                if pdf_links:
                    pdf_url = pdf_links[0].get_attribute('href')
                    st.write(f"✅ Found PDF URL: {pdf_url[:60]}...")
                    
                    # Download PDF
                    import requests
                    session = requests.Session()
                    
                    # Get cookies from browser
                    cookies = driver.get_cookies()
                    for cookie in cookies:
                        session.cookies.set(cookie['name'], cookie['value'])
                    
                    headers = {
                        'User-Agent': driver.execute_script("return navigator.userAgent;"),
                        'Referer': current_url
                    }
                    
                    response = session.get(pdf_url, headers=headers, timeout=30)
                    
                    if response.content.startswith(b'%PDF'):
                        st.success("✅ Successfully downloaded PDF!")
                        return response.content
                    else:
                        st.warning("⚠️ Downloaded content is not a PDF")
                else:
                    st.warning("⚠️ No PDF links found")
                        
            except Exception as e:
                st.error(f"❌ PDF download error: {str(e)}")
        
        return None
        
    except Exception as e:
        st.error(f"❌ Error with pre-authenticated session: {str(e)}")
        # Fall back to standard method
        st.warning("🔄 Falling back to standard authentication flow...")
        return download_doi_with_selenium(url)

def download_doi_with_selenium(url):
    """Download PDF from any DOI link using Selenium - Flexible approach"""
    driver = None
    try:
        st.write("🤖 Using Selenium for DOI PDF download (flexible mode)...")
        
        driver = create_selenium_driver()
        
        # Step 0: Establish institutional authentication session
        st.write(f"   0️⃣ Establishing institutional authentication...")
        try:
            # Visit ScienceDirect main page first to establish session
            driver.get("https://www.sciencedirect.com")
            time.sleep(3)
            
            # Check if we can access institutional content
            auth_test_url = "https://www.sciencedirect.com/science/article/pii/S0002944025001464"
            driver.get(auth_test_url)
            time.sleep(5)
            
            page_content = driver.page_source
            if 'Remote access' not in page_content and 'subscription' not in page_content.lower():
                st.write(f"   ✅ Institutional authentication established!")
            else:
                st.write(f"   ⚠️ May need institutional authentication")
                
        except Exception as auth_error:
            st.write(f"   ⚠️ Authentication setup error: {str(auth_error)[:50]}...")
        
        # Step 1: Visit DOI and handle any redirects
        st.write(f"   1️⃣ Following DOI redirect...")
        driver.get(url)
        
        # Wait for page load and any JavaScript redirects
        max_redirect_wait = 20
        start_time = time.time()
        previous_url = ""
        
        while time.time() - start_time < max_redirect_wait:
            time.sleep(2)
            current_url = driver.current_url
            
            if current_url == previous_url:
                # No more redirects
                break
            
            if current_url != url:
                st.write(f"   📍 Redirected to (full): {current_url}")
                if len(current_url) > 80:
                    st.write(f"   📍 Redirected to (short): {current_url[:80]}...")
            
            # Handle specific redirect scenarios
            if 'linkinghub.elsevier.com' in current_url:
                st.write(f"   ⏳ Waiting for Elsevier redirect...")
            elif 'onlinelibrary.wiley.com' in current_url:
                st.write(f"   ⏳ Waiting for Wiley redirect...")
            elif any(domain in current_url for domain in ['sciencedirect.com', 'springer.com', 'nature.com']):
                st.write(f"   ✅ Reached publisher site!")
                break
            
            previous_url = current_url
        
        # Get final URL
        final_url = driver.current_url
        st.write(f"   📍 Final destination (full): {final_url}")
        if len(final_url) > 80:
            st.write(f"   📍 Final destination (short): {final_url[:80]}...")
        
        # Step 2: Wait for dynamic content and analyze page structure
        st.write(f"   2️⃣ Waiting for dynamic content to load...")
        
        # Wait for JavaScript to load dynamic content - Extended for ScienceDirect
        max_wait_for_content = 30  # Increased from 15 to 30 seconds
        content_wait_start = time.time()
        
        while time.time() - content_wait_start < max_wait_for_content:
            try:
                # Check for PDF-related elements
                pdf_elements = driver.find_elements(By.XPATH, "//*[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'pdf')]")
                view_elements = driver.find_elements(By.XPATH, "//*[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'view')]")
                
                if pdf_elements or view_elements:
                    st.write(f"   ✅ Dynamic content loaded after {int(time.time() - content_wait_start)}s")
                    break
                    
                time.sleep(1)
                
            except Exception as wait_error:
                time.sleep(1)
                continue
        
        # Debug: Show page information after waiting
        try:
            page_title = driver.title
            st.write(f"   📄 Page title: {page_title[:80]}...")
            
            # Look for any buttons or links on the page
            all_buttons = driver.find_elements(By.TAG_NAME, 'button')
            all_links = driver.find_elements(By.TAG_NAME, 'a')
            st.write(f"   🔍 Found {len(all_buttons)} buttons and {len(all_links)} links on page")
            
            # Show some button texts for debugging
            if all_buttons:
                button_texts = []
                for btn in all_buttons[:10]:  # First 10 buttons
                    text = btn.text.strip()
                    if text and len(text) < 50:
                        button_texts.append(text)
                if button_texts:
                    st.write(f"   🔘 Button texts: {', '.join(button_texts[:5])}...")
            
            # Enhanced ScienceDirect debugging - Look for specific patterns
            if 'sciencedirect.com' in final_url:
                st.write(f"   🔬 ScienceDirect-specific analysis:")
                
                # Look for any elements containing "pdf", "view", "download" - enhanced for Download PDF
                pdf_related = driver.find_elements(By.XPATH, "//*[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'pdf') or contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'view') or contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'download')]")
                st.write(f"     📑 PDF/View/Download elements found: {len(pdf_related)}")
                
                # Specifically look for "Download PDF" links
                download_pdf_elements = driver.find_elements(By.XPATH, "//*[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'download pdf')]")
                if download_pdf_elements:
                    st.write(f"     🎯 'Download PDF' specific elements: {len(download_pdf_elements)}")
                    for i, elem in enumerate(download_pdf_elements[:3]):
                        try:
                            tag = elem.tag_name
                            text = elem.text.strip()[:30]
                            href = elem.get_attribute('href') if elem.tag_name == 'a' else 'N/A'
                            st.write(f"       {i+1}. {tag}: '{text}' href='{href[:40] if href != 'N/A' else 'N/A'}...'")
                        except:
                            continue
                
                # Show first few relevant elements
                for i, elem in enumerate(pdf_related[:5]):
                    try:
                        tag = elem.tag_name
                        text = elem.text.strip()[:30]
                        classes = elem.get_attribute('class') or ''
                        elem_id = elem.get_attribute('id') or ''
                        st.write(f"     {i+1}. {tag}: '{text}' class='{classes[:30]}' id='{elem_id[:20]}'")
                    except:
                        continue
                
                # Look for subscription/access related messages
                access_indicators = driver.find_elements(By.XPATH, "//*[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'access') or contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'subscription') or contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'sign in')]")
                if access_indicators:
                    st.write(f"     🔒 Access restriction indicators: {len(access_indicators)}")
                    for elem in access_indicators[:3]:
                        try:
                            text = elem.text.strip()[:50]
                            if text:
                                st.write(f"       ⚠️ '{text}'")
                        except:
                            continue
                
                # Check for alternative access options (Open Access, etc.)
                alt_access = driver.find_elements(By.XPATH, "//*[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'open access') or contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'free') or contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'available')]")
                if alt_access:
                    st.write(f"     🌍 Alternative access options: {len(alt_access)}")
                    for elem in alt_access[:3]:
                        try:
                            text = elem.text.strip()[:50]
                            if text and any(word in text.lower() for word in ['open', 'free', 'available']):
                                st.write(f"       ✅ '{text}'")
                        except:
                            continue
                
                # Check page source for hidden PDF patterns
                page_source = driver.page_source
                if 'pdfft' in page_source.lower():
                    st.write(f"     ✅ 'pdfft' pattern found in page source")
                if 'view pdf' in page_source.lower():
                    st.write(f"     ✅ 'view pdf' pattern found in page source")
                if 'download pdf' in page_source.lower():
                    st.write(f"     ✅ 'download pdf' pattern found in page source")
                
                # Detailed button/link analysis
                st.write(f"     📊 Detailed element analysis:")
                all_elements = driver.find_elements(By.XPATH, "//button | //a | //*[@role='button']")
                for i, elem in enumerate(all_elements[:15]):  # First 15 elements
                    try:
                        tag = elem.tag_name
                        text = elem.text.strip()[:30]
                        classes = elem.get_attribute('class') or ''
                        elem_id = elem.get_attribute('id') or ''
                        onclick = elem.get_attribute('onclick') or ''
                        data_aa = elem.get_attribute('data-aa-name') or ''
                        aria_label = elem.get_attribute('aria-label') or ''
                        href = elem.get_attribute('href') or ''
                        
                        # Show detailed info for potentially relevant elements
                        if any(keyword in (text + classes + elem_id + onclick + data_aa + aria_label + href).lower() 
                               for keyword in ['pdf', 'view', 'download', 'open', 'read']):
                            st.write(f"       🎯 {i+1}. {tag}: text='{text}' class='{classes[:20]}' id='{elem_id[:15]}' onclick='{onclick[:20]}' data-aa='{data_aa[:15]}' aria='{aria_label[:15]}' href='{href[:30]}'")
                        elif i < 10:  # Show first 10 elements regardless
                            st.write(f"       {i+1}. {tag}: '{text}' class='{classes[:20]}' id='{elem_id[:15]}'")
                    except Exception as elem_error:
                        continue
                
                # Deep DOM scan for ScienceDirect - Look for any clickable elements
                clickable_elements = driver.find_elements(By.XPATH, "//button | //a | //*[@onclick] | //*[@role='button'] | //*[@tabindex='0']")
                st.write(f"     🔍 Total clickable elements: {len(clickable_elements)}")
                
                # Check for elements with PDF-related attributes
                pdf_attributes = driver.find_elements(By.XPATH, "//*[contains(@data-aa-name, 'pdf') or contains(@data-testid, 'pdf') or contains(@aria-label, 'pdf') or contains(@title, 'pdf')]")
                if pdf_attributes:
                    st.write(f"     📎 PDF-related attributes found: {len(pdf_attributes)}")
                    for elem in pdf_attributes[:3]:
                        try:
                            tag = elem.tag_name
                            text = elem.text.strip()[:20]
                            data_aa = elem.get_attribute('data-aa-name') or ''
                            testid = elem.get_attribute('data-testid') or ''
                            aria_label = elem.get_attribute('aria-label') or ''
                            title = elem.get_attribute('title') or ''
                            st.write(f"       {tag}: '{text}' aa='{data_aa}' testid='{testid}' aria='{aria_label[:20]}' title='{title[:20]}'")
                        except:
                            continue
                
                # Wait and check for lazy-loaded content
                if time.time() - content_wait_start < 10:  # First 10 seconds, try scrolling
                    driver.execute_script("window.scrollTo(0, document.body.scrollHeight/3);")
                    time.sleep(1)
                    driver.execute_script("window.scrollTo(0, 0);")
                    st.write(f"     📜 Scrolled page to trigger lazy loading...")
                    
                    # Try clicking on potential interactive elements to trigger PDF loading
                    try:
                        # Force trigger any click events that might load PDF content
                        pdf_elements_count = driver.execute_script("""
                            // Comprehensive search for PDF-related elements
                            const elements = document.querySelectorAll('*');
                            let pdfElements = [];
                            let hiddenElements = [];
                            
                            elements.forEach(el => {
                                const text = el.textContent.toLowerCase();
                                const innerHTML = el.innerHTML.toLowerCase();
                                const attrs = Array.from(el.attributes).map(a => a.name + '=' + a.value).join(' ').toLowerCase();
                                
                                // Check for PDF-related content - enhanced for Download PDF
                                if (text.includes('pdf') || text.includes('view') || text.includes('download pdf') ||
                                    attrs.includes('pdf') || innerHTML.includes('pdf') || attrs.includes('view') || 
                                    attrs.includes('download') || text.includes('download')) {
                                    pdfElements.push({
                                        tag: el.tagName,
                                        text: el.textContent.trim().substring(0, 50),
                                        class: el.className,
                                        id: el.id,
                                        style: el.style.display,
                                        hidden: el.hidden,
                                        visible: el.offsetParent !== null
                                    });
                                }
                                
                                // Check for hidden elements that might contain View PDF or Download PDF
                                if (el.style.display === 'none' || el.hidden || el.offsetParent === null) {
                                    if (text.includes('pdf') || text.includes('view') || text.includes('download pdf') || 
                                        text.includes('download') || innerHTML.includes('pdf')) {
                                        hiddenElements.push({
                                            tag: el.tagName,
                                            text: el.textContent.trim().substring(0, 30),
                                            class: el.className,
                                            style: el.style.display
                                        });
                                    }
                                }
                            });
                            
                            // Try to trigger any pending JavaScript
                            window.dispatchEvent(new Event('load'));
                            window.dispatchEvent(new Event('DOMContentLoaded'));
                            
                            console.log('PDF elements found:', pdfElements);
                            console.log('Hidden PDF elements:', hiddenElements);
                            
                            return {
                                pdfElements: pdfElements,
                                hiddenElements: hiddenElements,
                                totalElements: elements.length
                            };
                        """)
                        
                        st.write(f"     🔍 Advanced JavaScript search completed...")
                        st.write(f"     📊 Total DOM elements: {pdf_elements_count.get('totalElements', 0)}")
                        
                        # Display found PDF elements
                        pdf_elements = pdf_elements_count.get('pdfElements', [])
                        if pdf_elements:
                            st.write(f"     🎯 PDF-related elements found: {len(pdf_elements)}")
                            for i, elem in enumerate(pdf_elements[:5]):
                                st.write(f"       {i+1}. {elem['tag']}: '{elem['text']}' class='{elem['class'][:20]}' visible={elem['visible']}")
                        
                        # Display hidden PDF elements
                        hidden_elements = pdf_elements_count.get('hiddenElements', [])
                        if hidden_elements:
                            st.write(f"     👻 Hidden PDF-related elements: {len(hidden_elements)}")
                            for i, elem in enumerate(hidden_elements[:3]):
                                st.write(f"       {i+1}. {elem['tag']}: '{elem['text']}' class='{elem['class'][:20]}' style='{elem['style']}'")
                                
                    except Exception as js_error:
                        st.write(f"     ⚠️ JavaScript search error: {str(js_error)[:50]}...")
                
                # Check for shadow DOM or iframe content
                try:
                    iframes = driver.find_elements(By.TAG_NAME, 'iframe')
                    if iframes:
                        st.write(f"     🖼️ Found {len(iframes)} iframes (potential PDF content)")
                        
                        # Inspect iframe content for PDF
                        for i, iframe in enumerate(iframes[:3]):  # Check first 3 iframes
                            try:
                                iframe_src = iframe.get_attribute('src')
                                iframe_id = iframe.get_attribute('id') or f'iframe_{i}'
                                st.write(f"       🔍 Iframe {i+1}: id='{iframe_id}' src='{iframe_src[:60] if iframe_src else 'None'}...'")
                                
                                # Check if iframe src contains PDF indicators
                                if iframe_src and any(pattern in iframe_src.lower() for pattern in ['pdf', 'pdfft', 'viewer']):
                                    st.write(f"       ✅ PDF-related iframe detected!")
                                    
                                    # Try to access iframe content
                                    try:
                                        driver.switch_to.frame(iframe)
                                        iframe_url = driver.current_url
                                        iframe_title = driver.title
                                        st.write(f"       📄 Iframe content: url='{iframe_url[:60]}...' title='{iframe_title[:30]}...'")
                                        
                                        # Look for PDF content inside iframe
                                        iframe_source = driver.page_source
                                        if 'application/pdf' in iframe_source or iframe_url.lower().endswith('.pdf'):
                                            st.write(f"       🎯 PDF content found in iframe!")
                                            # If this is a direct PDF URL, we could use it
                                            if iframe_url.lower().endswith('.pdf'):
                                                pdf_url = iframe_url
                                                pdf_found = True
                                                st.write(f"       ✅ Direct PDF URL: {pdf_url}")
                                        
                                        # Switch back to main frame
                                        driver.switch_to.default_content()
                                        
                                    except Exception as iframe_error:
                                        st.write(f"       ⚠️ Cannot access iframe content: {str(iframe_error)[:50]}...")
                                        driver.switch_to.default_content()  # Ensure we're back to main frame
                                        
                            except Exception as iframe_inspect_error:
                                st.write(f"       ❌ Iframe inspection error: {str(iframe_inspect_error)[:50]}...")
                                continue
                except:
                    pass
                        
        except Exception as debug_error:
            st.write(f"   ⚠️ Debug error: {str(debug_error)}")
        
        # Step 3: Flexible PDF search on the final page
        st.write(f"   3️⃣ Searching for PDF on the page...")
        
        # Strategy 1: Look for obvious PDF download links
        pdf_found = False
        pdf_url = None
        
        # Comprehensive PDF link search with ScienceDirect specific patterns
        pdf_search_patterns = [
            # Direct text searches - enhanced for Download PDF
            "//a[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'download pdf')]",
            "//button[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'download pdf')]",
            "//a[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'view pdf')]",
            "//button[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'view pdf')]",
            "//a[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'pdf')]",
            "//button[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'pdf')]",
            "//a[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'download')]",
            "//button[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'download')]",
            "//a[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'full text')]",
            
            # Attribute searches - enhanced for Download PDF
            "//a[contains(@href, 'pdf')]",
            "//a[contains(@href, '.pdf')]",
            "//a[contains(@href, 'pdfft')]",
            "//a[contains(@title, 'Download PDF')]",
            "//a[contains(@title, 'View PDF')]",
            "//a[contains(@title, 'PDF')]",
            "//a[contains(@aria-label, 'Download PDF')]",
            "//a[contains(@aria-label, 'View PDF')]",
            "//a[contains(@aria-label, 'PDF')]",
            "//button[contains(@title, 'Download PDF')]",
            "//button[contains(@title, 'View PDF')]",
            "//button[contains(@title, 'PDF')]",
            "//button[contains(@aria-label, 'Download PDF')]",
            "//button[contains(@aria-label, 'View PDF')]",
            "//button[contains(@aria-label, 'PDF')]",
            
            # Class and ID searches
            "//a[contains(@class, 'pdf')]",
            "//button[contains(@class, 'pdf')]",
            "//a[contains(@class, 'download')]",
            "//button[contains(@class, 'download')]",
            "//*[contains(@id, 'pdf')]",
            "//*[contains(@data-testid, 'pdf')]",
            
            # ScienceDirect specific patterns - enhanced for Download PDF
            "//button[contains(@class, 'ViewPDF')]",
            "//a[contains(@class, 'ViewPDF')]",
            "//button[contains(@class, 'DownloadPDF')]",
            "//a[contains(@class, 'DownloadPDF')]",
            "//button[contains(@data-aa-name, 'view-pdf')]",
            "//a[contains(@data-aa-name, 'view-pdf')]",
            "//button[contains(@data-aa-name, 'download-pdf')]",
            "//a[contains(@data-aa-name, 'download-pdf')]",
            "//*[@data-aa-name='view-pdf']",
            "//*[@data-aa-name='download-pdf']",
            "//*[@data-testid='view-pdf']",
            "//*[@data-testid='download-pdf']",
            "//button[contains(@aria-describedby, 'pdf')]",
            "//span[contains(text(), 'View PDF')]/parent::*/parent::a",
            "//span[contains(text(), 'View PDF')]/parent::*/parent::button",
            "//span[contains(text(), 'Download PDF')]/parent::*/parent::a",
            "//span[contains(text(), 'Download PDF')]/parent::*/parent::button",
            "//span[contains(text(), 'Download PDF')]/ancestor::a",
            "//span[contains(text(), 'Download PDF')]/ancestor::button",
            "//*[contains(@class, 'pdf-download')]",
            "//*[contains(@class, 'download-pdf')]",
            "//*[contains(@class, 'pdfLink')]"
        ]
        
        for i, pattern in enumerate(pdf_search_patterns):
            try:
                elements = driver.find_elements(By.XPATH, pattern)
                if elements:
                    st.write(f"   ✅ Found {len(elements)} elements with pattern {i+1}")
                    
                    for elem in elements[:3]:  # Check first 3 elements
                        href = elem.get_attribute('href')
                        text = elem.text.strip()
                        title = elem.get_attribute('title') or ''
                        
                        st.write(f"     📄 Element: '{text[:30]}...' href='{href[:50] if href else 'None'}...'")
                        
                        # If it's a link with href
                        if href and any(keyword in href.lower() for keyword in ['pdf', 'download', 'fulltext']):
                            pdf_url = href
                            st.write(f"   ✅ Found PDF link: {pdf_url[:80]}...")
                            pdf_found = True
                            break
                        
                        # If it's a clickable element (button), try clicking
                        elif elem.tag_name in ['button', 'a'] and any(keyword in text.lower() for keyword in ['pdf', 'download', 'full text']):
                            st.write(f"   🖱️ Clicking PDF element: '{text[:30]}...'")
                            try:
                                # Save current window handles
                                original_handles = driver.window_handles
                                elem.click()
                                time.sleep(3)
                                
                                # Check if new window/tab opened
                                new_handles = driver.window_handles
                                if len(new_handles) > len(original_handles):
                                    # New tab opened, switch to it
                                    driver.switch_to.window(new_handles[-1])
                                    new_url = driver.current_url
                                    st.write(f"   📝 New tab opened: {new_url[:80]}...")
                                    
                                    # Check if it's a PDF
                                    if new_url.lower().endswith('.pdf') or 'pdf' in new_url.lower():
                                        pdf_url = new_url
                                        pdf_found = True
                                        st.write(f"   ✅ PDF found in new tab!")
                                        break
                                    else:
                                        # Close new tab and return to original
                                        driver.close()
                                        driver.switch_to.window(original_handles[0])
                                else:
                                    # Check if current page changed to PDF
                                    current_url = driver.current_url
                                    if current_url != final_url:
                                        if current_url.lower().endswith('.pdf') or 'pdf' in current_url.lower():
                                            pdf_url = current_url
                                            pdf_found = True
                                            st.write(f"   ✅ Redirected to PDF!")
                                            break
                            except Exception as click_error:
                                st.write(f"     ⚠️ Click failed: {str(click_error)[:30]}...")
                                continue
                
                if pdf_found:
                    break
                    
            except Exception as pattern_error:
                continue
        
        # Strategy 2: If no direct PDF found, try alternative methods
        if not pdf_found:
            st.write(f"   4️⃣ Trying alternative PDF detection methods...")
            
            # Method 1: Enhanced "View PDF" and "Download PDF" detection for ScienceDirect
            try:
                # Multiple patterns for View PDF and Download PDF detection
                view_pdf_patterns = [
                    "//*[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'download pdf')]",
                    "//*[contains(translate(text(), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz'), 'view pdf')]",
                    "//button[contains(text(), 'Download PDF')]",
                    "//a[contains(text(), 'Download PDF')]",
                    "//button[contains(text(), 'View PDF')]",
                    "//a[contains(text(), 'View PDF')]", 
                    "//*[@aria-label='Download PDF']",
                    "//*[@aria-label='View PDF']",
                    "//*[@title='Download PDF']",
                    "//*[@title='View PDF']",
                    "//*[contains(@data-aa-name, 'download-pdf')]",
                    "//*[contains(@data-aa-name, 'view-pdf')]",
                    "//*[contains(@onclick, 'pdf')]",
                    "//button[contains(@class, 'pdf')]//span[contains(text(), 'Download')]",
                    "//a[contains(@class, 'pdf')]//span[contains(text(), 'Download')]",
                    "//button[contains(@class, 'pdf')]//span[contains(text(), 'View')]",
                    "//a[contains(@class, 'pdf')]//span[contains(text(), 'View')]"
                ]
                
                view_pdf_found = False
                for pattern in view_pdf_patterns:
                    if view_pdf_found:
                        break
                        
                    try:
                        view_pdf_elements = driver.find_elements(By.XPATH, pattern)
                        if view_pdf_elements:
                            st.write(f"   ✅ Found {len(view_pdf_elements)} PDF elements with pattern: {pattern[:50]}...")
                            
                            for elem in view_pdf_elements[:3]:
                                text = elem.text.strip()
                                tag = elem.tag_name
                                href = elem.get_attribute('href') if elem.tag_name == 'a' else 'N/A'
                                st.write(f"     📄 Element: {tag} '{text[:30]}...' href='{href[:40] if href != 'N/A' else 'N/A'}...'")
                                
                                try:
                                    # Save original window handles and tabs
                                    original_handles = driver.window_handles
                                    original_url = driver.current_url
                                    
                                    # Click the View PDF element
                                    st.write(f"     🖱️ Clicking View PDF element...")
                                    elem.click()
                                    time.sleep(5)  # Wait for response
                                    
                                    # Check for new window/tab
                                    new_handles = driver.window_handles
                                    if len(new_handles) > len(original_handles):
                                        st.write(f"     🆕 New tab opened!")
                                        driver.switch_to.window(new_handles[-1])
                                        new_url = driver.current_url
                                        st.write(f"     🔗 New tab URL: {new_url[:80]}...")
                                        
                                        # Check if it's the AWS S3 PDF URL
                                        if 'pdf.sciencedirectassets.com' in new_url or new_url.lower().endswith('.pdf'):
                                            pdf_url = new_url
                                            pdf_found = True
                                            view_pdf_found = True
                                            st.write(f"   ✅ ScienceDirect PDF URL found in new tab!")
                                            st.write(f"   🎯 URL: {pdf_url[:100]}...")
                                            break
                                        else:
                                            # Close new tab and return to original
                                            driver.close()
                                            driver.switch_to.window(original_handles[0])
                                    
                                    # Check if current page changed to PDF
                                    else:
                                        current_url = driver.current_url
                                        if current_url != original_url:
                                            st.write(f"     🔄 Page redirected: {current_url[:80]}...")
                                            if 'pdf.sciencedirectassets.com' in current_url or current_url.lower().endswith('.pdf'):
                                                pdf_url = current_url
                                                pdf_found = True
                                                view_pdf_found = True
                                                st.write(f"   ✅ ScienceDirect PDF URL found via redirect!")
                                                break
                                    
                                    # Check if iframe was populated
                                    try:
                                        iframes = driver.find_elements(By.TAG_NAME, 'iframe')
                                        for iframe in iframes:
                                            iframe_src = iframe.get_attribute('src')
                                            if iframe_src and ('pdf.sciencedirectassets.com' in iframe_src or '.pdf' in iframe_src):
                                                pdf_url = iframe_src
                                                pdf_found = True
                                                view_pdf_found = True
                                                st.write(f"   ✅ PDF URL found in iframe!")
                                                st.write(f"   🎯 URL: {pdf_url[:100]}...")
                                                break
                                    except:
                                        pass
                                        
                                except Exception as click_error:
                                    st.write(f"     ⚠️ Click failed: {str(click_error)[:50]}...")
                                    continue
                                    
                                if view_pdf_found:
                                    break
                    except:
                        continue
                        
            except Exception as view_pdf_error:
                st.write(f"   ⚠️ View PDF detection error: {str(view_pdf_error)[:50]}...")
            
            # Method 2: Look for any downloadable content
            if not pdf_found:
                st.write(f"   5️⃣ Searching for any downloadable content...")
                
                # Look for any links that might lead to PDF
                all_links = driver.find_elements(By.TAG_NAME, 'a')
                potential_pdf_links = []
                
                for link in all_links:
                    href = link.get_attribute('href')
                    text = link.text.strip().lower()
                    
                    if href and any(keyword in (href.lower() + ' ' + text) for keyword in 
                        ['full', 'view', 'read', 'access', 'open', 'article', 'content', 'text']):
                        potential_pdf_links.append((link, href, text))
                
                if potential_pdf_links:
                    st.write(f"   📋 Found {len(potential_pdf_links)} potential content links")
                    
                    # Try the most promising links
                    for link, href, text in potential_pdf_links[:3]:
                        st.write(f"   🔍 Trying: '{text[:30]}...' -> {href[:50]}...")
                        try:
                            original_url = driver.current_url
                            link.click()
                            time.sleep(3)
                            
                            new_url = driver.current_url
                            if new_url.lower().endswith('.pdf') or 'application/pdf' in driver.page_source:
                                pdf_url = new_url
                                pdf_found = True
                                st.write(f"   ✅ Found PDF content!")
                                break
                            else:
                                # Go back if not PDF
                                driver.back()
                                time.sleep(2)
                        except:
                            continue
        
        # Detect publisher for specific handling
        publisher = 'unknown'
        if 'wiley.com' in final_url:
            publisher = 'wiley'
        elif 'sciencedirect.com' in final_url:
            publisher = 'elsevier'
        elif 'springer.com' in final_url:
            publisher = 'springer'
        elif 'nature.com' in final_url:
            publisher = 'nature'
        elif 'mdpi.com' in final_url:
            publisher = 'mdpi'
        elif 'pmc.ncbi.nlm.nih.gov' in final_url:
            publisher = 'pmc'
        elif 'atsjournals.org' in final_url:
            publisher = 'ats'
        
        st.write(f"   🏢 Publisher detected: {publisher}")
        
        # Fallback: Try direct PDF URL construction for ScienceDirect
        if not pdf_found and 'sciencedirect.com' in final_url:
            st.write(f"   6️⃣ Fallback: Trying direct PDF URL construction...")
            
            # Extract PII from URL
            import re
            pii_match = re.search(r'/pii/([A-Z0-9]+)', final_url)
            if pii_match:
                pii = pii_match.group(1)
                
                # Try multiple PDF URL patterns
                pdf_patterns = [
                    f"https://www.sciencedirect.com/science/article/pii/{pii}/pdfft",
                    f"https://www.sciencedirect.com/science/article/pii/{pii}/pdf",
                    f"https://pdf.sciencedirectassets.com/main.pdf?_ob=ArticleURL&_method=download&_eid=1-s2.0-{pii}&origin=article",
                    f"https://api.elsevier.com/content/article/pii/{pii}?view=FULL"
                ]
                
                # Test each pattern
                for i, test_url in enumerate(pdf_patterns):
                    st.write(f"     {i+1}. Testing (full URL): {test_url}")
                    
                    try:
                        # Navigate to PDF URL directly
                        driver.get(test_url)
                        time.sleep(3)
                        
                        # For /pdfft URLs, extract AWS S3 URL directly from JavaScript
                        if '/pdfft' in test_url:
                            st.write(f"       🔍 Extracting AWS S3 URL from ScienceDirect JavaScript...")
                            
                            # Wait for JavaScript to execute
                            st.write(f"       ⏳ Waiting for JavaScript redirect...")
                            time.sleep(3)
                            
                            # Check if page has redirected
                            current_url = driver.current_url
                            if 'pdf.sciencedirectassets.com' in current_url:
                                st.write(f"       ✅ Page redirected to AWS S3!")
                                test_url = current_url
                            else:
                                # Extract PDF URL from page source
                                try:
                                    page_source = driver.page_source
                                    
                                    # Execute JavaScript to get the redirect URL
                                    try:
                                        js_result = driver.execute_script("""
                                            // Look for window.location assignments in scripts
                                            var scripts = document.getElementsByTagName('script');
                                            for (var i = 0; i < scripts.length; i++) {
                                                var content = scripts[i].innerHTML;
                                                var match = content.match(/window\\.location\\s*=\\s*["']([^"']*pdf\\.sciencedirectassets\\.com[^"']*)/);
                                                if (match) {
                                                    return match[1];
                                                }
                                            }
                                            
                                            // Check noscript content
                                            var noscripts = document.getElementsByTagName('noscript');
                                            for (var i = 0; i < noscripts.length; i++) {
                                                var content = noscripts[i].innerHTML;
                                                var match = content.match(/href=["']([^"']*pdf\\.sciencedirectassets\\.com[^"']*)/);
                                                if (match) {
                                                    return match[1];
                                                }
                                            }
                                            
                                            return null;
                                        """)
                                        
                                        if js_result:
                                            st.write(f"       ✅ Found URL via JavaScript extraction")
                                            extracted_url = js_result
                                        else:
                                            st.write(f"       ⚠️ JavaScript extraction returned null")
                                            extracted_url = None
                                    except:
                                        extracted_url = None
                                    
                                    # Fallback to regex patterns
                                    if not extracted_url:
                                        import re
                                        
                                        # Pattern 1: window.location assignment (with various quote styles)
                                        location_patterns = [
                                            r"window\.location\s*=\s*['\"]([^'\"]*pdf\.sciencedirectassets\.com[^'\"]*)['\"]",
                                            r"window\.location\s*=\s*`([^`]*pdf\.sciencedirectassets\.com[^`]*)`",
                                            r"window\.location\.href\s*=\s*['\"]([^'\"]*pdf\.sciencedirectassets\.com[^'\"]*)['\"]"
                                        ]
                                        
                                        for pattern in location_patterns:
                                            match = re.search(pattern, page_source)
                                            if match:
                                                extracted_url = match.group(1)
                                                st.write(f"       ✅ Found URL in window.location pattern")
                                                break
                                        
                                        # Pattern 2: noscript href
                                        if not extracted_url:
                                            noscript_pattern = r'<a\s+href=["\']([^"\']*pdf\.sciencedirectassets\.com[^"\']*)["\']'
                                            noscript_match = re.search(noscript_pattern, page_source)
                                            if noscript_match:
                                                extracted_url = noscript_match.group(1)
                                                st.write(f"       ✅ Found URL in noscript link")
                                        
                                        # Pattern 3: Any AWS S3 URL in the page
                                        if not extracted_url:
                                            general_pattern = r'(https?://pdf\.sciencedirectassets\.com/[^"\'\s<>]*)'
                                            general_match = re.search(general_pattern, page_source)
                                            if general_match:
                                                extracted_url = general_match.group(1)
                                                st.write(f"       ✅ Found AWS S3 URL in page source")
                                    
                                    if extracted_url:
                                        # Clean up any HTML entities
                                        import html
                                        extracted_url = html.unescape(extracted_url)
                                        
                                        # Ensure it's a complete URL
                                        if not extracted_url.startswith('http'):
                                            extracted_url = 'https:' + extracted_url
                                        
                                        st.write(f"       🎯 EXTRACTED PDF URL: {extracted_url[:100]}...")
                                        st.write(f"       📍 Full URL: {extracted_url}")
                                        
                                        # Navigate to the extracted URL
                                        driver.get(extracted_url)
                                        time.sleep(3)
                                        test_url = extracted_url  # Update for further processing
                                    else:
                                        st.write(f"       ❌ No AWS S3 URL found in page source")
                                    
                                    # Debug: Show page source snippets
                                    st.write(f"       🔍 Debug: Analyzing page content...")
                                    
                                    # Check for "window.location" patterns
                                    window_location_snippets = []
                                    for match in re.finditer(r'window\.location[^;]{0,200}', page_source):
                                        snippet = match.group(0)
                                        window_location_snippets.append(snippet[:100])
                                    
                                    if window_location_snippets:
                                        st.write(f"       📋 Found window.location patterns:")
                                        for i, snippet in enumerate(window_location_snippets[:3]):
                                            st.write(f"         {i+1}. {snippet}...")
                                    
                                    # Check for any ScienceDirect asset URLs
                                    asset_urls = re.findall(r'sciencedirectassets\.com[^"\'\s]*', page_source)
                                    if asset_urls:
                                        st.write(f"       📋 Found ScienceDirect asset references:")
                                        for i, url in enumerate(asset_urls[:3]):
                                            st.write(f"         {i+1}. ...{url[:80]}...")
                                    
                                    # Check for AWS references
                                    aws_refs = re.findall(r'[^"\'\s]*X-Amz[^"\'\s]*', page_source)
                                    if aws_refs:
                                        st.write(f"       📋 Found AWS S3 signature references:")
                                        for i, ref in enumerate(aws_refs[:2]):
                                            st.write(f"         {i+1}. {ref[:60]}...")
                                    
                                    # Check page structure
                                    if 'Are you a robot' in page_source:
                                        st.write(f"       🤖 DETECTED: 'Are you a robot?' challenge page")
                                    if 'redirectUser' in page_source:
                                        st.write(f"       🔄 DETECTED: redirectUser function in page")
                                    if 'cloudflare' in page_source.lower():
                                        st.write(f"       ☁️ DETECTED: Cloudflare protection")
                                    
                                    # Show first 500 characters of page for debugging
                                    st.write(f"       📄 Page start: {page_source[:500]}...")
                                    
                                except Exception as extract_error:
                                    st.write(f"       ⚠️ URL extraction failed: {str(extract_error)[:50]}...")
                                
                                # Fallback: wait for automatic redirect
                                st.write(f"       ⏳ Falling back to automatic redirect detection...")
                                for redirect_wait in range(8):  # Reduced to 16 seconds
                                    time.sleep(2)
                                    current_check = driver.current_url
                                    
                                    if 'pdf.sciencedirectassets.com' in current_check:
                                        st.write(f"       🎯 REDIRECT SUCCESS! New URL: {current_check}")
                                        test_url = current_check
                                        break
                                    
                                    if redirect_wait % 3 == 2:
                                        st.write(f"       ⏳ Still waiting... ({(redirect_wait+1)*2}s)")
                        
                        current_content_type = None
                        try:
                            # Check content type via JavaScript
                            current_content_type = driver.execute_script("return document.contentType;")
                        except:
                            pass
                        
                        current_url_after = driver.current_url
                        st.write(f"       Current URL (full): {current_url_after}")
                        st.write(f"       Content-Type: {current_content_type}")
                        
                        # Check if we got PDF content
                        if current_content_type and 'pdf' in current_content_type.lower():
                            st.write(f"       ✅ PDF content type detected!")
                            pdf_url = current_url_after
                            pdf_found = True
                            break
                        elif current_url_after.lower().endswith('.pdf'):
                            st.write(f"       ✅ PDF URL format detected!")
                            pdf_url = current_url_after
                            pdf_found = True
                            break
                        elif 'pdf.sciencedirectassets.com' in current_url_after:
                            st.write(f"       ✅ ScienceDirect assets detected!")
                            pdf_url = current_url_after
                            pdf_found = True
                            
                            # Try to get PDF content directly from current page
                            try:
                                current_page_source = driver.page_source
                                if current_page_source.startswith('%PDF') or len(current_page_source) > 100000:
                                    st.write(f"       📄 PDF content available in current page!")
                                    
                                    # Page source might be the PDF content
                                    if current_page_source.startswith('%PDF'):
                                        st.success("   ✅ PDF content extracted from page source!")
                                        return current_page_source.encode('latin-1')
                                    else:
                                        st.write(f"       📊 Page size: {len(current_page_source)} chars - checking for binary content...")
                                        
                                        # Check if large page contains PDF-like content
                                        if len(current_page_source) > 500000:  # Large content, likely PDF
                                            st.write(f"       🔍 Large content detected - attempting PDF extraction...")
                                            
                                            # Look for PDF magic number or structure in the content
                                            try:
                                                # Try to encode as binary and check for PDF header with Unicode handling
                                                try:
                                                    potential_pdf = current_page_source.encode('latin-1')
                                                except UnicodeEncodeError:
                                                    # Handle Unicode characters
                                                    st.write(f"       🔧 Unicode detected, cleaning content...")
                                                    clean_content = current_page_source.replace('\u200b', '').replace('\ufeff', '')
                                                    try:
                                                        potential_pdf = clean_content.encode('latin-1')
                                                    except UnicodeEncodeError:
                                                        # Raw character conversion as fallback
                                                        potential_pdf = bytes([ord(c) & 0xFF for c in current_page_source])
                                                
                                                # Check various positions for PDF header (sometimes there's HTML wrapper)
                                                pdf_positions = []
                                                search_content = potential_pdf[:50000]  # Search first 50KB
                                                
                                                for i in range(0, len(search_content) - 4, 100):
                                                    if search_content[i:i+4] == b'%PDF':
                                                        pdf_positions.append(i)
                                                
                                                if pdf_positions:
                                                    st.write(f"       ✅ Found PDF header at position(s): {pdf_positions}")
                                                    # Extract from first PDF position
                                                    pdf_start = pdf_positions[0]
                                                    pdf_content = potential_pdf[pdf_start:]
                                                    
                                                    if len(pdf_content) > 10000:  # Reasonable PDF size
                                                        st.success("   ✅ PDF content extracted from large page source!")
                                                        return pdf_content
                                                
                                                # Fallback: check if content has PDF-like structure
                                                if (b'/Type' in potential_pdf[:10000] and 
                                                    b'/Catalog' in potential_pdf[:10000] and 
                                                    len(potential_pdf) > 100000):
                                                    st.write(f"       🎯 PDF structure detected without magic number")
                                                    st.success("   ✅ PDF-like content extracted!")
                                                    return potential_pdf
                                                    
                                            except Exception as binary_error:
                                                st.write(f"       ⚠️ Binary extraction failed: {str(binary_error)[:50]}...")
                                        
                                        # Simplified JavaScript approach for smaller content
                                        if len(current_page_source) < 500000:
                                            try:
                                                # Simpler JavaScript without complex fetch
                                                pdf_binary = driver.execute_script("""
                                                    try {
                                                        var content = document.documentElement.innerHTML;
                                                        if (content.length > 100000) {
                                                            return content;
                                                        }
                                                        return null;
                                                    } catch(e) {
                                                        return null;
                                                    }
                                                """)
                                                
                                                if pdf_binary and len(pdf_binary) > 100000:
                                                    try:
                                                        potential_pdf = pdf_binary.encode('latin-1')
                                                        if b'%PDF' in potential_pdf[:1000] or b'/Type' in potential_pdf[:5000]:
                                                            st.success("   ✅ PDF content extracted via simplified JavaScript!")
                                                            return potential_pdf
                                                    except:
                                                        pass
                                            except Exception as simple_js_error:
                                                st.write(f"       ⚠️ Simplified JavaScript failed: {str(simple_js_error)[:50]}...")
                            except Exception as content_error:
                                st.write(f"       ⚠️ Content extraction error: {str(content_error)[:50]}...")
                            
                            break
                        else:
                            # Check page source for PDF indicators
                            page_source = driver.page_source
                            if 'PDF' in page_source[:500] and len(page_source) < 2000:
                                st.write(f"       🔍 Possible PDF content (small page with PDF text)")
                                pdf_url = current_url_after
                                pdf_found = True
                                break
                            else:
                                st.write(f"       ❌ Not PDF content (Page size: {len(page_source)} chars)")
                        
                    except Exception as pattern_error:
                        st.write(f"       ⚠️ Pattern {i+1} failed: {str(pattern_error)[:50]}...")
                        continue
                
                if pdf_found:
                    st.write(f"   ✅ Direct PDF URL construction successful!")
                else:
                    st.write(f"   ❌ All PDF URL patterns failed")
            else:
                st.write(f"   ❌ Could not extract PII from URL")
        
        # Step 3: Download the PDF
        if pdf_found and pdf_url:
            st.write(f"   4️⃣ Downloading PDF from (full URL): {pdf_url}")
        if len(pdf_url) > 80:
            st.write(f"   4️⃣ Downloading PDF from (short): {pdf_url[:80]}...")
            
            # Get cookies and user agent from browser session
            cookies = driver.get_cookies()
            user_agent = driver.execute_script("return navigator.userAgent;")
            
            # Create session with browser cookies
            session = requests.Session()
            for cookie in cookies:
                session.cookies.set(
                    cookie['name'], 
                    cookie['value'],
                    domain=cookie.get('domain'),
                    path=cookie.get('path', '/')
                )
            
            # Prepare headers that mimic the browser
            headers = {
                'User-Agent': user_agent,
                'Accept': 'application/pdf,text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8',
                'Accept-Language': 'en-US,en;q=0.9',
                'Accept-Encoding': 'gzip, deflate, br',
                'DNT': '1',
                'Connection': 'keep-alive',
                'Upgrade-Insecure-Requests': '1',
                'Sec-Fetch-Dest': 'document',
                'Sec-Fetch-Mode': 'navigate',
                'Sec-Fetch-Site': 'same-origin',
                'Cache-Control': 'max-age=0',
                'Referer': final_url
            }
            
            # Enhanced download with multiple attempts
            download_success = False
            download_attempts = [
                # Attempt 1: Use session with all cookies and referrer
                {
                    'method': 'session_with_referrer',
                    'headers': {**headers, 'Referer': final_url},
                    'stream': False
                },
                # Attempt 2: Direct browser navigation and content extraction
                {
                    'method': 'browser_navigation',
                    'headers': headers,
                    'stream': False
                },
                # Attempt 3: Session with modified headers
                {
                    'method': 'session_modified',
                    'headers': {
                        **headers,
                        'Accept': 'application/pdf,*/*',
                        'Cache-Control': 'no-cache',
                        'Pragma': 'no-cache'
                    },
                    'stream': True
                }
            ]
            
            for attempt_num, attempt in enumerate(download_attempts, 1):
                if download_success:
                    break
                    
                st.write(f"   📥 Download attempt {attempt_num}: {attempt['method']}")
                
                try:
                    if attempt['method'] == 'browser_navigation':
                        # Try navigating to PDF URL with browser and extracting content
                        st.write(f"     🌐 Navigating browser to PDF URL...")
                        driver.get(pdf_url)
                        time.sleep(5)
                        
                        # Check if browser loaded PDF directly
                        current_content_type = None
                        try:
                            current_content_type = driver.execute_script("return document.contentType;")
                        except:
                            pass
                            
                        # Check if PDF is displayed in browser (regardless of content-type)
                        current_url_check = driver.current_url
                        page_source = driver.page_source
                        
                        st.write(f"     📊 Browser state (full URL): {current_url_check}")
                        st.write(f"     📊 Content-Type: {current_content_type}, Page size: {len(page_source)} chars")
                        
                        # Check for PDF viewer in browser
                        if (current_content_type and 'pdf' in current_content_type.lower()) or \
                           'pdf.sciencedirectassets.com' in current_url_check or \
                           len(page_source) > 500000:  # Large content suggests PDF
                            
                            st.write(f"     ✅ PDF detected in browser!")
                            
                            # Method 1: Try direct page source extraction (for embedded PDFs)
                            if page_source.startswith('%PDF'):
                                st.success("     🎯 PDF found directly in page source!")
                                return page_source.encode('latin-1')
                            
                            # Method 2: Check if large content contains PDF
                            elif len(page_source) > 500000:
                                st.write(f"     🔍 Large content detected - checking for PDF...")
                                try:
                                    # Enhanced encoding for large content
                                    try:
                                        potential_pdf = page_source.encode('latin-1')
                                    except UnicodeEncodeError:
                                        st.write(f"     🔧 Unicode in large content, using raw conversion...")
                                        potential_pdf = bytes([ord(c) & 0xFF for c in page_source])
                                    
                                    # Look for PDF magic number
                                    if b'%PDF' in potential_pdf[:1000]:
                                        st.success("     ✅ PDF header found in large content!")
                                        return potential_pdf
                                    elif b'/Type' in potential_pdf[:5000] and b'/Catalog' in potential_pdf[:5000]:
                                        st.success("     ✅ PDF structure detected!")
                                        return potential_pdf
                                except Exception as extract_error:
                                    st.write(f"     ⚠️ Large content extraction failed: {str(extract_error)[:50]}...")
                            
                            # Method 3: JavaScript-based extraction with timeout
                            try:
                                st.write(f"     🔧 Attempting JavaScript extraction...")
                                
                                # Simple and robust JavaScript approach
                                pdf_content = driver.execute_script("""
                                    // Check if we're on a PDF page
                                    var contentType = document.contentType || '';
                                    var url = window.location.href;
                                    
                                    if (contentType.toLowerCase().includes('pdf') || 
                                        url.includes('pdf.sciencedirectassets.com')) {
                                        
                                        // Try to access the PDF content
                                        try {
                                            var xhr = new XMLHttpRequest();
                                            xhr.open('GET', url, false);  // Synchronous request
                                            xhr.overrideMimeType('text/plain; charset=x-user-defined');
                                            xhr.send();
                                            
                                            if (xhr.status === 200) {
                                                return xhr.responseText;
                                            }
                                        } catch(e) {
                                            console.log('XHR failed:', e);
                                        }
                                        
                                        // Fallback: return page content if large
                                        var content = document.documentElement.outerHTML || document.body.innerHTML;
                                        if (content.length > 100000) {
                                            return content;
                                        }
                                    }
                                    
                                    return null;
                                """)
                                
                                if pdf_content and len(pdf_content) > 10000:
                                    st.write(f"     📦 JavaScript returned {len(pdf_content)} characters")
                                    
                                    # Enhanced encoding handling for Unicode characters
                                    encoding_attempts = [
                                        ('latin-1', 'Direct binary encoding'),
                                        ('utf-8', 'UTF-8 with binary conversion'),
                                        ('utf-16', 'UTF-16 encoding'),
                                        ('raw', 'Raw character processing')
                                    ]
                                    
                                    for encoding_name, encoding_desc in encoding_attempts:
                                        try:
                                            st.write(f"       🔧 Trying {encoding_desc}...")
                                            
                                            if encoding_name == 'raw':
                                                # Raw character processing - convert each character to byte
                                                pdf_bytes = bytes([ord(c) & 0xFF for c in pdf_content])
                                            elif encoding_name == 'utf-8':
                                                # UTF-8 with special handling
                                                # Remove problematic Unicode characters first
                                                clean_content = pdf_content.replace('\u200b', '').replace('\ufeff', '')
                                                pdf_bytes = clean_content.encode('latin-1')
                                            elif encoding_name == 'utf-16':
                                                # Try UTF-16 and convert to bytes
                                                utf16_bytes = pdf_content.encode('utf-16le')
                                                # Take every other byte (simple conversion)
                                                pdf_bytes = bytes(utf16_bytes[i] for i in range(0, len(utf16_bytes), 2))
                                            else:
                                                # Standard latin-1
                                                pdf_bytes = pdf_content.encode(encoding_name)
                                            
                                            # Check if the result looks like PDF
                                            if pdf_bytes.startswith(b'%PDF'):
                                                st.success(f"     ✅ PDF content extracted via {encoding_desc}!")
                                                return pdf_bytes
                                            elif len(pdf_bytes) > 100000 and (b'/Type' in pdf_bytes[:5000] or b'PDF' in pdf_bytes[:1000]):
                                                st.success(f"     ✅ PDF structure detected via {encoding_desc}!")
                                                return pdf_bytes
                                            else:
                                                st.write(f"       ❌ {encoding_desc} - doesn't produce valid PDF")
                                                # Show first few bytes for debugging
                                                first_bytes = pdf_bytes[:20] if len(pdf_bytes) >= 20 else pdf_bytes
                                                st.write(f"       📊 First bytes: {first_bytes}")
                                                
                                        except Exception as encoding_error:
                                            st.write(f"       ⚠️ {encoding_desc} failed: {str(encoding_error)[:50]}...")
                                            continue
                                    
                                    # If all encoding attempts failed, try content analysis
                                    st.write(f"     🔍 Analyzing content structure...")
                                    
                                    # Look for PDF patterns in the string content itself
                                    if '%PDF' in pdf_content[:1000]:
                                        st.write(f"       ✅ Found %PDF in content!")
                                        try:
                                            # Extract from the %PDF position
                                            pdf_start = pdf_content.find('%PDF')
                                            pdf_part = pdf_content[pdf_start:]
                                            # Try raw conversion on just the PDF part
                                            pdf_bytes = bytes([ord(c) & 0xFF for c in pdf_part])
                                            if pdf_bytes.startswith(b'%PDF'):
                                                st.success("     ✅ PDF extracted from content substring!")
                                                return pdf_bytes
                                        except Exception as substring_error:
                                            st.write(f"       ⚠️ Substring extraction failed: {str(substring_error)[:50]}...")
                                    
                                    st.write(f"     ❌ All encoding methods failed to produce valid PDF")
                                else:
                                    st.write(f"     ❌ JavaScript extraction returned insufficient content")
                                    
                            except Exception as js_error:
                                st.write(f"     ⚠️ JavaScript extraction failed: {str(js_error)[:50]}...")
                        else:
                            st.write(f"     ❌ No PDF content detected (Content-Type: {current_content_type})")
                    
                    else:
                        # Regular session download
                        response = session.get(
                            pdf_url, 
                            headers=attempt['headers'], 
                            timeout=30, 
                            allow_redirects=True,
                            stream=attempt['stream']
                        )
                        
                        st.write(f"     📊 Response: {response.status_code}, Content-Type: {response.headers.get('Content-Type', 'Unknown')}")
                        
                        if response.status_code == 200:
                            content = response.content
                            if content.startswith(b'%PDF'):
                                st.success(f"   ✅ Successfully downloaded PDF with attempt {attempt_num}!")
                                download_success = True
                                return content
                            else:
                                st.write(f"     ⚠️ Content is not PDF (Size: {len(content)} bytes)")
                                # Check if content might be PDF despite headers
                                if len(content) > 10000 and (b'PDF' in content[:1000] or b'/Type /Catalog' in content[:5000]):
                                    st.info(f"     💡 Content appears to be PDF despite headers")
                                    download_success = True
                                    return content
                        elif response.status_code == 403:
                            st.write(f"     🔒 Access denied (403) - trying next method...")
                        else:
                            st.write(f"     ❌ HTTP Error {response.status_code}")
                            
                except Exception as download_error:
                    st.write(f"     ❌ Attempt {attempt_num} failed: {str(download_error)[:50]}...")
                    continue
            
            if not download_success:
                st.error("   ❌ All download attempts failed")
        
        else:
            st.warning("   ⚠️ No PDF download link found on the page")
            
            # Fallback: Try to save current page if it appears to be PDF content
            current_content_type = driver.execute_script("return document.contentType;")
            if 'pdf' in current_content_type.lower():
                st.info("   💡 Current page appears to be PDF content")
                try:
                    # Get page source as potential PDF
                    page_source = driver.page_source
                    if page_source.startswith('%PDF'):
                        st.success("   ✅ Captured PDF from current page!")
                        return page_source.encode()
                except:
                    pass
            
    except Exception as e:
        st.error(f"   ❌ Flexible Selenium error: {str(e)}")
    finally:
        if driver:
            driver.quit()
    
    return None


def find_pdf_by_selectors(driver):
    """Find PDF links using common selectors"""
    pdf_selectors = [
        "a[href*='pdf']",
        "a[href*='PDF']",
        "a[title*='PDF']",
        "a[aria-label*='PDF']",
        "button[aria-label*='PDF']",
        "a.pdf-download",
        "a.pdf-link",
        ".pdf-download-btn",
        "a[data-pdf-url]",
        "button.download-pdf",
        "[class*='pdf-download']",
        "[class*='download-pdf']"
    ]
    
    for selector in pdf_selectors:
        try:
            elements = driver.find_elements(By.CSS_SELECTOR, selector)
            for elem in elements:
                href = elem.get_attribute('href')
                if href and ('pdf' in href.lower() or '/doi/pdf/' in href):
                    return href
        except:
            continue
    
    return None


def find_pdf_in_all_links(driver):
    """Search all links for PDF"""
    links = driver.find_elements(By.TAG_NAME, 'a')
    for link in links:
        try:
            href = link.get_attribute('href')
            text = link.text.lower()
            if href and ('pdf' in href.lower() or 'pdf' in text or 'download' in text):
                return href
        except:
            continue
    
    return None


def click_download_button(driver):
    """Try to click download button and get PDF URL"""
    download_selectors = [
        "button:contains('Download')",
        "a:contains('Download')",
        "[aria-label*='Download']",
        ".download-button"
    ]
    
    for selector in download_selectors:
        try:
            if ':contains(' in selector:
                xpath = f"//*[contains(text(), 'Download')]"
                elem = driver.find_element(By.XPATH, xpath)
            else:
                elem = driver.find_element(By.CSS_SELECTOR, selector)
            
            elem.click()
            time.sleep(2)
            
            # Check if new window/tab opened
            if len(driver.window_handles) > 1:
                driver.switch_to.window(driver.window_handles[-1])
                return driver.current_url
                
        except:
            continue
    
    return None


def construct_pdf_url(driver, current_url, publisher):
    """Construct PDF URL based on publisher patterns"""
    if publisher == 'wiley':
        # Wiley pattern: /doi/10.xxx/xxx -> /doi/pdf/10.xxx/xxx
        if '/doi/' in current_url:
            return current_url.replace('/doi/', '/doi/pdf/')
    
    elif publisher == 'elsevier':
        # Extract PII and construct PDF URL
        import re
        pii_match = re.search(r'/pii/([A-Z0-9]+)', current_url)
        if pii_match:
            pii = pii_match.group(1)
            return f"https://www.sciencedirect.com/science/article/pii/{pii}/pdfft"
    
    elif publisher == 'springer':
        # Springer pattern
        if '/article/' in current_url:
            return current_url.replace('/article/', '/content/pdf/') + '.pdf'
    
    elif publisher == 'nature':
        # Nature pattern
        if '/articles/' in current_url:
            return current_url + '.pdf'
    
    return None


def download_with_selenium(url, site_type='auto'):
    """
    Generic Selenium download function that detects site type
    
    Args:
        url: The URL to download from
        site_type: 'auto', 'pmc', 'sciencedirect', 'doi', etc.
    """
    # Always use Selenium for DOI links
    if 'doi.org' in url:
        return download_doi_with_selenium(url)
    
    if site_type == 'auto':
        # Auto-detect site type
        if 'pmc.ncbi.nlm.nih.gov' in url:
            site_type = 'pmc'
        elif 'sciencedirect.com' in url:
            site_type = 'sciencedirect'
        elif 'elsevier.com' in url:
            site_type = 'sciencedirect'
        else:
            # For any unknown site, try generic DOI approach
            return download_doi_with_selenium(url)
    
    if site_type == 'pmc':
        return download_pmc_with_selenium(url)
    elif site_type == 'sciencedirect':
        return download_sciencedirect_with_selenium(url)
    elif site_type == 'doi':
        return download_doi_with_selenium(url)
    else:
        # Default to DOI approach
        return download_doi_with_selenium(url)