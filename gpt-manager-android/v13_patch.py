from pathlib import Path

root=Path('.')
java=root/'app/src/main/java/com/gptyar/manager/MainActivity.java'
index=root/'app/src/main/assets/index.html'
build=root/'app/build.gradle'
s=java.read_text(encoding='utf-8')

s=s.replace('import android.webkit.CookieManager;','import android.webkit.CookieManager;\nimport android.webkit.WebChromeClient;')
s=s.replace('GPTManagerAndroid/1.2','GPTManagerAndroid/1.3')

# New bridge actions.
old='''                    case "adminNotificationDelete": return adminNotificationDelete(body.optLong("notificationId", -1)).toString();
                    default: return error("درخواست ناشناخته است.").toString();'''
new='''                    case "adminNotificationDelete": return adminNotificationDelete(body.optLong("notificationId", -1)).toString();
                    case "resume": return resumeSession().toString();
                    case "deactivate": return deactivate().toString();
                    case "openActiveService": return openActiveService().toString();
                    default: return error("درخواست ناشناخته است.").toString();'''
if old not in s: raise SystemExit('v13 bridge marker missing')
s=s.replace(old,new)

# Persist token history synchronously so it is visible immediately after logout.
s=s.replace('.edit().putString("recent_tokens_v1", next.toString()).apply();','.edit().putString("recent_tokens_v1", next.toString()).commit();')
s=s.replace('.edit().remove("recent_tokens_v1").apply();','.edit().remove("recent_tokens_v1").commit();')

# Make WebView act like an embedded browser and never finish Activity on window.close.
old='''        s.setDomStorageEnabled(true);
        s.setDatabaseEnabled(true);'''
new='''        s.setDomStorageEnabled(true);
        s.setDatabaseEnabled(true);
        s.setSupportMultipleWindows(false);
        s.setJavaScriptCanOpenWindowsAutomatically(true);'''
if old not in s: raise SystemExit('v13 web settings marker missing')
s=s.replace(old,new)

marker='''        webView.setWebViewClient(new WebViewClient() {'''
chrome='''        webView.setWebChromeClient(new WebChromeClient() {
            @Override public void onCloseWindow(WebView window) {
                // Sites sometimes call window.close() after auth/popups. Keep the app alive.
                returnToDashboard();
            }
        });

''' + marker
if marker not in s: raise SystemExit('v13 WebViewClient marker missing')
s=s.replace(marker,chrome)

# More natural back behavior: navigate browser history first; dashboard back backgrounds app instead of destroying it.
start=s.index('    @Override public void onBackPressed() {')
end=s.index('    private String deviceId()', start)
newback='''    @Override public void onBackPressed() {
        String u = webView.getUrl();
        if (u != null && !u.startsWith("file:///android_asset/")) {
            if (webView.canGoBack()) webView.goBack(); else returnToDashboard();
            return;
        }
        moveTaskToBack(true);
    }

'''
s=s[:start]+newback+s[end:]

# Add active-state helper and session resume.
marker='''    private JSONObject logout() throws Exception {'''
helpers='''    private JSONObject activeState() {
        SharedPreferences p=getSharedPreferences("gpt_manager", MODE_PRIVATE);
        try {
            return new JSONObject()
                    .put("subscriptionId", p.getLong("active_subscription_id", 0L))
                    .put("agentType", p.getString("active_agent_type", ""))
                    .put("target", p.getString("active_target", ""));
        } catch (JSONException e) { return new JSONObject(); }
    }

    private void saveActiveState(long subscriptionId, String agentType, String target) {
        getSharedPreferences("gpt_manager", MODE_PRIVATE).edit()
                .putLong("active_subscription_id", subscriptionId)
                .putString("active_agent_type", agentType == null ? "" : agentType)
                .putString("active_target", target == null ? "" : target)
                .commit();
    }

    private void clearActiveState() {
        getSharedPreferences("gpt_manager", MODE_PRIVATE).edit()
                .remove("active_subscription_id").remove("active_agent_type").remove("active_target").commit();
    }

    private JSONObject resumeSession() throws Exception {
        String token, role;
        synchronized (stateLock) { token=sessionToken; role=sessionRole; }
        if (token.isEmpty() || role.isEmpty()) return error("نشست فعالی وجود ندارد.");
        if ("admin".equals(role)) {
            JSONObject subs=adminSubscribers();
            return ok().put("role","admin").put("token",token).put("subscribers",subs.opt("data"))
                    .put("profile",new JSONObject()).put("deviceId",deviceId()).put("active",activeState());
        }
        JSONObject o=userOverview();
        if (!o.optBoolean("ok",false)) return o;
        return o.put("role","user").put("token",token).put("deviceId",deviceId()).put("active",activeState());
    }

    private JSONObject deactivate() throws Exception {
        clearServiceState();
        clearActiveState();
        return ok();
    }

    private JSONObject openActiveService() throws Exception {
        JSONObject a=activeState();
        String target=a.optString("target","");
        if (target.isEmpty()) return error("حساب فعالی وجود ندارد.");
        main.post(() -> { webView.removeJavascriptInterface("AndroidApi"); webView.loadUrl(target); });
        return ok();
    }

''' + marker
if marker not in s: raise SystemExit('v13 logout marker missing')
s=s.replace(marker,helpers)

# logout clears active browser session but keeps recent-token history.
old='''    private JSONObject logout() throws Exception {
        clearServiceState();
        synchronized (stateLock) { sessionToken = ""; sessionRole = ""; cachedSubscribers = new JSONArray(); }
        return ok();
    }'''
new='''    private JSONObject logout() throws Exception {
        clearServiceState();
        clearActiveState();
        synchronized (stateLock) { sessionToken = ""; sessionRole = ""; cachedSubscribers = new JSONArray(); }
        return ok();
    }'''
if old not in s: raise SystemExit('v13 logout body marker missing')
s=s.replace(old,new)

# Include active state in login/overview results.
s=s.replace('return ok().put("role", "admin").put("profile", r.data).put("subscribers", subs.opt("data")).put("deviceId", deviceId());',
            'return ok().put("role", "admin").put("profile", r.data).put("subscribers", subs.opt("data")).put("deviceId", deviceId()).put("active", activeState());')
s=s.replace('overview.put("role", "user").put("deviceId", deviceId());',
            'overview.put("role", "user").put("deviceId", deviceId()).put("active", activeState());')

# When activation succeeds, persist which subscription is active just like extension activeAgents/activeSubscription state.
old='''                String target = serviceUrl(ref.type);
                main.post(() -> {
                    webView.removeJavascriptInterface("AndroidApi");
                    webView.loadUrl(target);
                });
                return ok().put("activated", true).put("agentType", ref.type).put("subscriptionId", subscriptionId).put("diagnostic", diag);'''
new='''                String target = serviceUrl(ref.type);
                saveActiveState(subscriptionId, ref.type, target);
                main.post(() -> {
                    webView.removeJavascriptInterface("AndroidApi");
                    webView.loadUrl(target);
                });
                return ok().put("activated", true).put("agentType", ref.type).put("target", target).put("subscriptionId", subscriptionId).put("diagnostic", diag);'''
if old not in s: raise SystemExit('v13 activation success marker missing')
s=s.replace(old,new)

# If activation fully fails, preserve no stale active UI.
old='''        clearServiceState();
        return activationError(diag, last == null ? "فعال‌سازی ناموفق بود." : (last.getMessage()==null?last.toString():last.getMessage()));'''
new='''        clearServiceState();
        clearActiveState();
        return activationError(diag, last == null ? "فعال‌سازی ناموفق بود." : (last.getMessage()==null?last.toString():last.getMessage()));'''
if old not in s: raise SystemExit('v13 activation failure marker missing')
s=s.replace(old,new)

required=['new WebChromeClient()','onCloseWindow(WebView window)','case "resume"','case "deactivate"','active_subscription_id','saveActiveState(subscriptionId, ref.type, target)','commit();']
for x in required:
    if x not in s: raise SystemExit('v13 missing '+x)
java.write_text(s,encoding='utf-8')
index.write_text((root/'v13_index.html').read_text(encoding='utf-8'),encoding='utf-8')

b=build.read_text(encoding='utf-8')
b=b.replace("versionCode 3\n        versionName '1.2.0'","versionCode 4\n        versionName '1.3.0'")
if "versionName '1.3.0'" not in b: raise SystemExit('v13 version update failed')
build.write_text(b,encoding='utf-8')
print('Applied GPT Manager Android v1.3.0 final patch')
