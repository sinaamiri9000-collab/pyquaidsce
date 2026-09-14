from pathlib import Path

root=Path('.')
java=root/'app/src/main/java/com/gptyar/manager/MainActivity.java'
index=root/'app/src/main/assets/index.html'
build=root/'app/build.gradle'

s=java.read_text(encoding='utf-8')
s=s.replace('import android.os.Looper;','import android.os.Looper;\nimport android.os.Message;')
s=s.replace('import android.widget.Button;','import android.widget.Button;\nimport android.widget.FrameLayout;')
s=s.replace('import java.util.Base64;','import java.util.ArrayDeque;\nimport java.util.Base64;')
s=s.replace('GPTManagerAndroid/1.3','GPTManagerAndroid/1.4')

# Browser container and popup stack. Closing a site popup must never finish the Activity.
old='''    private WebView webView;\n    private LinearLayout toolbar;\n    private TextView toolbarTitle;'''
new='''    private WebView webView;\n    private FrameLayout webContainer;\n    private final ArrayDeque<WebView> popupViews = new ArrayDeque<>();\n    private LinearLayout toolbar;\n    private TextView toolbarTitle;'''
if old not in s: raise SystemExit('v14 fields marker missing')
s=s.replace(old,new)

# Multiple windows are handled by child WebViews instead of letting window.close affect the Activity.
s=s.replace('s.setSupportMultipleWindows(false);','s.setSupportMultipleWindows(true);')

old='''        webView.setWebChromeClient(new WebChromeClient() {\n            @Override public void onCloseWindow(WebView window) {\n                // Sites sometimes call window.close() after auth/popups. Keep the app alive.\n                returnToDashboard();\n            }\n        });'''
new='''        webView.setWebChromeClient(new WebChromeClient() {\n            @Override public boolean onCreateWindow(WebView view, boolean isDialog, boolean isUserGesture, Message resultMsg) {\n                return createPopupWebView(resultMsg);\n            }\n            @Override public void onCloseWindow(WebView window) {\n                if (window != null && window != webView) closePopupWebView(window);\n                else returnToDashboard();\n            }\n        });'''
if old not in s: raise SystemExit('v14 chrome client marker missing')
s=s.replace(old,new)

old='''        root.addView(toolbar, new LinearLayout.LayoutParams(LinearLayout.LayoutParams.MATCH_PARENT, LinearLayout.LayoutParams.WRAP_CONTENT));\n        root.addView(webView, new LinearLayout.LayoutParams(LinearLayout.LayoutParams.MATCH_PARENT, 0, 1));\n        setContentView(root);'''
new='''        root.addView(toolbar, new LinearLayout.LayoutParams(LinearLayout.LayoutParams.MATCH_PARENT, LinearLayout.LayoutParams.WRAP_CONTENT));\n        webContainer = new FrameLayout(this);\n        webContainer.addView(webView, new FrameLayout.LayoutParams(FrameLayout.LayoutParams.MATCH_PARENT, FrameLayout.LayoutParams.MATCH_PARENT));\n        root.addView(webContainer, new LinearLayout.LayoutParams(LinearLayout.LayoutParams.MATCH_PARENT, 0, 1));\n        setContentView(root);'''
if old not in s: raise SystemExit('v14 root webview marker missing')
s=s.replace(old,new)

# Replace simple dashboard/back behavior with an embedded-browser lifecycle.
start=s.index('    private String hostOf(String url) {')
end=s.index('    private String deviceId()', start)
block='''    private String hostOf(String url) {\n        try { return new URL(url).getHost(); } catch (Exception e) { return "سرویس"; }\n    }\n\n    @SuppressLint("SetJavaScriptEnabled")\n    private void configurePopupWebView(WebView child) {\n        WebSettings ps = child.getSettings();\n        ps.setJavaScriptEnabled(true);\n        ps.setDomStorageEnabled(true);\n        ps.setDatabaseEnabled(true);\n        ps.setAllowFileAccess(false);\n        ps.setAllowContentAccess(false);\n        ps.setMixedContentMode(WebSettings.MIXED_CONTENT_NEVER_ALLOW);\n        ps.setSupportMultipleWindows(true);\n        ps.setJavaScriptCanOpenWindowsAutomatically(true);\n        ps.setUserAgentString(ps.getUserAgentString() + " GPTManagerAndroid/1.4");\n        CookieManager.getInstance().setAcceptCookie(true);\n        CookieManager.getInstance().setAcceptThirdPartyCookies(child, true);\n        child.setWebViewClient(new WebViewClient() {\n            @Override public boolean shouldOverrideUrlLoading(WebView view, WebResourceRequest request) { return false; }\n            @Override public void onPageFinished(WebView view, String url) {\n                toolbar.setVisibility(View.VISIBLE);\n                toolbarTitle.setText(hostOf(url));\n            }\n            @Override public void onReceivedHttpAuthRequest(WebView view, HttpAuthHandler handler, String host, String realm) {\n                String u, p;\n                synchronized (stateLock) { u=currentProxyUser; p=currentProxyPass; }\n                if (!u.isEmpty()) handler.proceed(u,p); else handler.cancel();\n            }\n            @Override public void onReceivedSslError(WebView view, SslErrorHandler handler, SslError error) { handler.cancel(); }\n        });\n        child.setWebChromeClient(new WebChromeClient() {\n            @Override public boolean onCreateWindow(WebView view, boolean isDialog, boolean isUserGesture, Message resultMsg) {\n                return createPopupWebView(resultMsg);\n            }\n            @Override public void onCloseWindow(WebView window) { closePopupWebView(window); }\n        });\n    }\n\n    private boolean createPopupWebView(Message resultMsg) {\n        if (resultMsg == null || webContainer == null) return false;\n        try {\n            WebView child = new WebView(this);\n            configurePopupWebView(child);\n            popupViews.addLast(child);\n            webContainer.addView(child, new FrameLayout.LayoutParams(FrameLayout.LayoutParams.MATCH_PARENT, FrameLayout.LayoutParams.MATCH_PARENT));\n            child.bringToFront();\n            toolbar.setVisibility(View.VISIBLE);\n            WebView.WebViewTransport transport = (WebView.WebViewTransport) resultMsg.obj;\n            transport.setWebView(child);\n            resultMsg.sendToTarget();\n            return true;\n        } catch (Throwable t) {\n            return false;\n        }\n    }\n\n    private void closePopupWebView(WebView child) {\n        if (child == null || child == webView) { returnToDashboard(); return; }\n        popupViews.remove(child);\n        try { webContainer.removeView(child); } catch (Throwable ignored) {}\n        try { child.stopLoading(); child.loadUrl("about:blank"); child.destroy(); } catch (Throwable ignored) {}\n        WebView top = popupViews.peekLast();\n        if (top != null) {\n            top.bringToFront();\n            toolbar.setVisibility(View.VISIBLE);\n            toolbarTitle.setText(hostOf(top.getUrl()));\n        } else {\n            webView.bringToFront();\n            String u=webView.getUrl();\n            boolean local=u!=null && u.startsWith("file:///android_asset/");\n            toolbar.setVisibility(local?View.GONE:View.VISIBLE);\n            toolbarTitle.setText(local?"":hostOf(u));\n        }\n    }\n\n    private void closeAllPopupWebViews() {\n        while (!popupViews.isEmpty()) {\n            WebView child=popupViews.pollLast();\n            try { webContainer.removeView(child); } catch (Throwable ignored) {}\n            try { child.stopLoading(); child.loadUrl("about:blank"); child.destroy(); } catch (Throwable ignored) {}\n        }\n        if (webView != null) webView.bringToFront();\n    }\n\n    private WebView topWebView() {\n        WebView top=popupViews.peekLast();\n        return top==null?webView:top;\n    }\n\n    private void returnToDashboard() {\n        closeAllPopupWebViews();\n        webView.addJavascriptInterface(nativeBridge, "AndroidApi");\n        webView.loadUrl(DASHBOARD);\n    }\n\n    private void browserBack() {\n        WebView top=topWebView();\n        if (top != null && top != webView) {\n            if (top.canGoBack()) top.goBack(); else closePopupWebView(top);\n            return;\n        }\n        String u=webView==null?null:webView.getUrl();\n        if (u != null && !u.startsWith("file:///android_asset/")) {\n            if (webView.canGoBack()) webView.goBack(); else returnToDashboard();\n            return;\n        }\n        moveTaskToBack(true);\n    }\n\n    @Override public void onBackPressed() { browserBack(); }\n\n    @Override protected void onDestroy() {\n        closeAllPopupWebViews();\n        clearProxy();\n        try { proxyExecutor.shutdownNow(); } catch (Throwable ignored) {}\n        super.onDestroy();\n    }\n\n'''
s=s[:start]+block+s[end:]

# Logout returns the already-persisted token history in the same bridge response.
old='''    private JSONObject logout() throws Exception {\n        clearServiceState();\n        clearActiveState();\n        synchronized (stateLock) { sessionToken = ""; sessionRole = ""; cachedSubscribers = new JSONArray(); }\n        return ok();\n    }'''
new='''    private JSONObject logout() throws Exception {\n        clearServiceState();\n        clearActiveState();\n        synchronized (stateLock) { sessionToken = ""; sessionRole = ""; cachedSubscribers = new JSONArray(); }\n        return ok().put("recentTokens", recentTokens());\n    }'''
if old not in s: raise SystemExit('v14 logout marker missing')
s=s.replace(old,new)

# Keep UI state synchronized on ordinary refreshes too.
old='''.put("subscriptionLoadOk", agents.ok)\n                .put("subscriptionError", agents.ok ? JSONObject.NULL : messageOf(agents));'''
new='''.put("subscriptionLoadOk", agents.ok)\n                .put("subscriptionError", agents.ok ? JSONObject.NULL : messageOf(agents))\n                .put("active", activeState());'''
if old not in s: raise SystemExit('v14 overview marker missing')
s=s.replace(old,new)

# Match the extension's switch semantics: activating a new subscription starts from a clean previous state.
marker='''        step(steps, "subscription", true, ref.type + " / " + ref.accountName);\n\n        JSONObject pb = new JSONObject().put("token", token);'''
replacement='''        step(steps, "subscription", true, ref.type + " / " + ref.accountName);\n        clearServiceState();\n        clearActiveState();\n\n        JSONObject pb = new JSONObject().put("token", token);'''
if marker not in s: raise SystemExit('v14 activation clean marker missing')
s=s.replace(marker,replacement)

# HTML/UI fixes.
h=index.read_text(encoding='utf-8')
h=h.replace('v1.3.0','v1.4.0')
h=h.replace('.subscription{margin-top:10px;padding:13px;border-radius:15px;background:#0e2639;border:1px solid #1c455c}', '.subscription{margin-top:10px;padding:13px;border-radius:15px;background:#0e2639;border:1px solid #1c455c}.subscription.active-sub{border-color:#1c7d61;box-shadow:0 0 0 1px #1bc58d22 inset;background:linear-gradient(180deg,#0f2b37,#0d2433)}')

old="""function renderRecent(){const r=call('recentTokens');recentTokens=r.ok?(r.data||[]):[];const w=$('#recentWrap'),list=$('#recentList');if(!recentTokens.length){w.classList.add('hidden');list.innerHTML='';return}w.classList.remove('hidden');list.innerHTML=recentTokens.map((x,i)=>`<div class=\"recent\"><div class=\"recentMain\" onclick=\"loginRecent(${i})\"><div class=\"recentLabel\">${esc(x.label||((x.role==='admin')?'مدیر':'کاربر'))}</div><div class=\"recentToken\">${esc(maskToken(x.token))}</div><div class=\"recentTime\">آخرین ورود: ${pdate(Number(x.lastUsed||0),true)}</div></div><button class=\"iconBtn\" onclick=\"forgetRecent(event,${i})\">×</button></div>`).join('')}"""
new="""function renderRecent(data=null){if(Array.isArray(data))recentTokens=data;else{const r=call('recentTokens');recentTokens=r.ok?(r.data||[]):[]}const w=$('#recentWrap'),list=$('#recentList');if(!recentTokens.length){w.classList.add('hidden');list.innerHTML='';return}w.classList.remove('hidden');list.innerHTML=recentTokens.map((x,i)=>`<div class=\"recent\"><div class=\"recentMain\" onclick=\"loginRecent(${i})\"><div class=\"recentLabel\">${esc(x.label||((x.role==='admin')?'مدیر':'کاربر'))}</div><div class=\"recentToken\">${esc(maskToken(x.token))}</div><div class=\"recentTime\">آخرین ورود: ${pdate(Number(x.lastUsed||0),true)}</div></div><button class=\"iconBtn\" onclick=\"forgetRecent(event,${i})\">×</button></div>`).join('')}"""
if old not in h: raise SystemExit('v14 renderRecent marker missing')
h=h.replace(old,new)

old="""function logout(){const previous=state.token;call('logout');state={role:null,data:null,token:null,activeSubscriptionId:0,activeAgentType:'',activeTarget:''};adminRows=[];adminNotices=[];$('#appView').classList.add('hidden');$('#loginView').classList.remove('hidden');$('#tokenInput').value='';content.innerHTML='';renderRecent();setTimeout(renderRecent,80)}"""
new="""function logout(){const r=call('logout');state={role:null,data:null,token:null,activeSubscriptionId:0,activeAgentType:'',activeTarget:''};adminRows=[];adminNotices=[];$('#appView').classList.add('hidden');$('#loginView').classList.remove('hidden');$('#tokenInput').value='';content.innerHTML='';const rows=r&&r.ok&&Array.isArray(r.recentTokens)?r.recentTokens:null;renderRecent(rows);setTimeout(()=>renderRecent(),120)}"""
if old not in h: raise SystemExit('v14 JS logout marker missing')
h=h.replace(old,new)

old="""function refreshCurrent(){const b=$('#refreshBtn');busy(b,true);setTimeout(()=>{if(state.role==='admin'){const r=call('adminSubscribers');if(r.ok){adminRows=r.data||[];if(adminTab==='notices')loadAdminNotices(false);renderAdmin();toast('اطلاعات بروزرسانی شد')}else toast(r.error,true)}else{const r=call('overview');if(r.ok){r.role='user';r.deviceId=state.data?.deviceId;state.data=r;renderUser(r);toast('اشتراک‌ها بروزرسانی شدند')}else toast(r.error,true)}busy(b,false)},20)}"""
new="""function refreshCurrent(){const b=$('#refreshBtn');busy(b,true);setTimeout(()=>{if(state.role==='admin'){const r=call('adminSubscribers');if(r.ok){adminRows=r.data||[];if(adminTab==='notices')loadAdminNotices(false);renderAdmin();toast('اطلاعات بروزرسانی شد')}else toast(r.error,true)}else{const r=call('overview');if(r.ok){r.role='user';r.deviceId=state.data?.deviceId;state.data=r;const a=r.active||{};state.activeSubscriptionId=Number(a.subscriptionId||state.activeSubscriptionId||0);state.activeAgentType=a.agentType||state.activeAgentType||'';state.activeTarget=a.target||state.activeTarget||'';renderUser(r);toast('اشتراک‌ها بروزرسانی شدند')}else toast(r.error,true)}busy(b,false)},20)}"""
if old not in h: raise SystemExit('v14 refresh marker missing')
h=h.replace(old,new)

start=h.index('function serviceHtml(a){')
end=h.index('function activate(',start)
service="""function serviceHtml(a){const type=agentType(a),subs=subscriptions(a);return `<div class=\"service\"><div class=\"serviceHead\"><div class=\"serviceName\">${esc(a.title||a.name||type)}</div><span class=\"badge\">${esc(type)}</span></div>${subs.length?subs.map(s=>{const id=subId(s),active=Number(state.activeSubscriptionId)===Number(id);return `<div class=\"subscription ${active?'active-sub':''}\" id=\"sub-${id}\"><div class=\"subTop\"><div><div class=\"subName\">${esc(acct(s))}</div><div class=\"subMeta\">شناسه ${fa(id)} • انقضا ${pdate(s.expiresAt||s.expirationDate)}</div></div><span class=\"badge\" style=\"${active?'color:#43d9a8;border-color:#21745d':'color:#8fa8b9'}\">${active?'فعال':'غیرفعال'}</span></div><div class=\"actions\">${active?`<button class=\"btn green small\" onclick=\"openActiveService(this)\">ورود به سرویس</button><button class=\"btn red small\" onclick=\"deactivateCurrent(this)\">غیرفعال‌سازی</button>`:`<button class=\"btn primary small\" onclick=\"activate(${id},this)\">فعال‌سازی و ورود</button>`}</div><div class=\"diag hidden\" data-diag></div></div>`}).join(''):'<div class=\"empty\">اشتراکی برای این سرویس ثبت نشده است.</div>'}</div>`}\n"""
h=h[:start]+service+h[end:]
h=h.replace("state.activeTarget=r.target||'';if(state.data)renderUser(state.data);toast('حساب فعال شد')", "state.activeTarget=r.target||'';if(state.data)renderUser(state.data);toast('اشتراک فعال شد')")
h=h.replace("if(!autoResume())renderRecent();", "if(!autoResume()){renderRecent();setTimeout(()=>renderRecent(),120)}")

for x in ['function renderRecent(data=null)','recentTokens)?r.recentTokens:null','فعال‌سازی و ورود','غیرفعال‌سازی','v1.4.0']:
    if x not in h: raise SystemExit('v14 HTML missing '+x)
index.write_text(h,encoding='utf-8')

b=build.read_text(encoding='utf-8')
b=b.replace("versionCode 4\n        versionName '1.3.0'","versionCode 5\n        versionName '1.4.0'")
if "versionName '1.4.0'" not in b: raise SystemExit('v14 version update failed')
build.write_text(b,encoding='utf-8')

required=['createPopupWebView(Message resultMsg)','closePopupWebView(WebView child)','popupViews','return ok().put("recentTokens", recentTokens())','.put("active", activeState())','clearServiceState();\n        clearActiveState();']
for x in required:
    if x not in s: raise SystemExit('v14 Java missing '+x)
java.write_text(s,encoding='utf-8')
print('Applied GPT Manager Android v1.4.0 final fixes')
