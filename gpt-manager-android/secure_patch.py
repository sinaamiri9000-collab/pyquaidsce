from pathlib import Path

p = Path('app/src/main/java/com/gptyar/manager/MainActivity.java')
s = p.read_text(encoding='utf-8')

s = s.replace(
    '    private JSONArray cachedSubscribers = new JSONArray();\n',
    '    private JSONArray cachedSubscribers = new JSONArray();\n    private final NativeBridge nativeBridge = new NativeBridge();\n'
)
s = s.replace(
    '        webView.addJavascriptInterface(new NativeBridge(), "AndroidApi");',
    '        webView.addJavascriptInterface(nativeBridge, "AndroidApi");'
)
s = s.replace(
    '            @Override public boolean shouldOverrideUrlLoading(WebView view, WebResourceRequest request) { return false; }',
    '''            @Override public boolean shouldOverrideUrlLoading(WebView view, WebResourceRequest request) {\n                String next = request == null || request.getUrl() == null ? "" : request.getUrl().toString();\n                if (!next.startsWith("file:///android_asset/")) view.removeJavascriptInterface("AndroidApi");\n                return false;\n            }'''
)
s = s.replace(
    '    private void returnToDashboard() { webView.loadUrl(DASHBOARD); }',
    '''    private void returnToDashboard() {\n        webView.addJavascriptInterface(nativeBridge, "AndroidApi");\n        webView.loadUrl(DASHBOARD);\n    }'''
)
s = s.replace(
    '''                String target = serviceUrl(ref.type);\n                main.post(() -> webView.loadUrl(target));''',
    '''                String target = serviceUrl(ref.type);\n                main.post(() -> {\n                    webView.removeJavascriptInterface("AndroidApi");\n                    webView.loadUrl(target);\n                });'''
)

required = [
    'private final NativeBridge nativeBridge = new NativeBridge();',
    'view.removeJavascriptInterface("AndroidApi")',
    'webView.removeJavascriptInterface("AndroidApi")',
]
for marker in required:
    if marker not in s:
        raise SystemExit(f'security patch missing marker: {marker}')

p.write_text(s, encoding='utf-8')
print('Applied WebView bridge isolation patch')
