from pathlib import Path

root = Path('.')
java = root / 'app/src/main/java/com/gptyar/manager/MainActivity.java'
index = root / 'app/src/main/assets/index.html'
manifest = root / 'app/src/main/AndroidManifest.xml'
build = root / 'app/build.gradle'
icon = root / 'app/src/main/res/drawable/ic_launcher.xml'

s = java.read_text(encoding='utf-8')
s = s.replace('GPTManagerAndroid/1.0', 'GPTManagerAndroid/1.1')

old = '''                    case "adminCreate": return adminCreate(body).toString();
                    default: return error("درخواست ناشناخته است.").toString();'''
new = '''                    case "adminCreate": return adminCreate(body).toString();
                    case "recentTokens": return recentTokensResult().toString();
                    case "forgetRecentToken": return forgetRecentToken(body.optString("token")).toString();
                    case "clearRecentTokens": return clearRecentTokens().toString();
                    default: return error("درخواست ناشناخته است.").toString();'''
if old not in s: raise SystemExit('v11 patch: bridge switch marker missing')
s = s.replace(old, new)

start = s.index('    private JSONObject login(String token) throws Exception {')
end = s.index('    private JSONObject logout() throws Exception {')
login = '''    private JSONObject login(String token) throws Exception {
        token = token == null ? "" : token.trim();
        if (token.length() < 4) return error("توکن معتبر وارد کنید.");

        clearProxy();
        synchronized (stateLock) {
            sessionToken = ""; sessionRole = ""; cachedSubscribers = new JSONArray();
            currentProxyUser = ""; currentProxyPass = "";
        }

        if (token.startsWith("@")) {
            JSONObject b = new JSONObject().put("token", token).put("version", VERSION);
            ApiResponse r = api("/accupdator/validate", "POST", token, b, null, true);
            if (!r.ok || r.data.optString("message").isEmpty() || !r.data.has("account") || !r.data.has("agent"))
                return upstreamError(r, "توکن مدیریتی معتبر نیست.");
            synchronized (stateLock) { sessionToken = token; sessionRole = "admin"; }
            JSONObject subs = adminSubscribers();
            JSONObject account = r.data.optJSONObject("account");
            String label = account == null ? "" : firstString(account, "name", "username", "title");
            saveRecentToken(token, "admin", label.isEmpty() ? "مدیر" : label);
            return ok().put("role", "admin").put("profile", r.data).put("subscribers", subs.opt("data")).put("deviceId", deviceId());
        }

        JSONObject b = new JSONObject()
                .put("token", token).put("device_id", deviceId()).put("version", VERSION)
                .put("client", new JSONObject().put("product", "PLUGIN").put("version", VERSION).put("platform", "CHROME"));
        ApiResponse v = api("/client/validate", "POST", token, b, null, true);
        if (!v.ok && (v.status == 408 || v.status == 429 || v.status >= 500)) {
            Thread.sleep(650);
            v = api("/client/validate", "POST", token, b, null, true);
        }
        JSONObject notice = v.data.optJSONObject("notice");
        if (notice != null) return error(notice.optString("message", notice.optString("text", "ورود توسط سرور متوقف شد."))).put("notice", notice);
        if (!v.ok || !v.data.optBoolean("valid", false)) return upstreamError(v, "توکن معتبر نیست.");

        synchronized (stateLock) { sessionToken = token; sessionRole = "user"; }
        JSONObject overview = userOverview();
        if (!overview.optBoolean("ok", false)) {
            synchronized (stateLock) { sessionToken = ""; sessionRole = ""; }
            return overview;
        }
        overview.put("role", "user").put("deviceId", deviceId());
        JSONObject user = overview.optJSONObject("user");
        String label = user == null ? "" : firstString(user, "name", "username", "email", "title");
        saveRecentToken(token, "user", label.isEmpty() ? "کاربر" : label);
        return overview;
    }

'''
s = s[:start] + login + s[end:]

old = '''    private JSONObject userOverview() throws Exception {
        String token;
        synchronized (stateLock) { token = sessionToken; }
        if (token.isEmpty()) return error("ابتدا وارد شوید.");
        ApiResponse info = api("/clients/user-info", "POST", token, null, null, true);
        ApiResponse agents = api("/clients/subscription/agents-with-subscriptions", "GET", token, null, null, true);
        ApiResponse notes = api("/clients/notifications/active", "GET", token, null, null, true);
        return ok().put("user", info.ok ? info.data : JSONObject.NULL)
                .put("agents", agents.ok ? normalizeArray(agents.data, "agents", "data") : new JSONArray())
                .put("notifications", notes.ok ? normalizeArray(notes.data, "data") : new JSONArray());
    }
'''
new = '''    private JSONObject userOverview() throws Exception {
        String token;
        synchronized (stateLock) { token = sessionToken; }
        if (token.isEmpty()) return error("ابتدا وارد شوید.");
        ApiResponse info = api("/clients/user-info", "POST", token, null, null, true);
        ApiResponse agents = api("/clients/subscription/agents-with-subscriptions", "GET", token, null, null, true);
        ApiResponse notes = api("/clients/notifications/active", "GET", token, null, null, true);

        JSONArray agentRows = agents.ok ? normalizeArray(agents.data, "agents", "data", "result") : new JSONArray();
        JSONArray noteRows = notes.ok ? normalizeArray(notes.data, "notifications", "data", "result") : new JSONArray();
        if (!info.ok && !agents.ok) {
            return error("اطلاعات حساب دریافت نشد. دوباره تلاش کنید.")
                    .put("userError", messageOf(info)).put("subscriptionsError", messageOf(agents));
        }
        return ok().put("user", info.ok ? info.data : JSONObject.NULL)
                .put("agents", agentRows)
                .put("notifications", noteRows)
                .put("subscriptionLoadOk", agents.ok)
                .put("subscriptionError", agents.ok ? JSONObject.NULL : messageOf(agents));
    }
'''
if old not in s: raise SystemExit('v11 patch: overview marker missing')
s = s.replace(old, new)

old = '''        JSONObject data;
        try { data = text.isEmpty() ? new JSONObject() : new JSONObject(text); }
        catch (JSONException ex) { data = new JSONObject().put("message", text); }'''
new = '''        JSONObject data;
        String trimmed = text == null ? "" : text.trim();
        try {
            if (trimmed.isEmpty()) data = new JSONObject();
            else if (trimmed.startsWith("[")) data = new JSONObject().put("_root", new JSONArray(trimmed));
            else data = new JSONObject(trimmed);
        } catch (JSONException ex) { data = new JSONObject().put("message", text); }'''
if old not in s: raise SystemExit('v11 patch: API parser marker missing')
s = s.replace(old, new)

old = '''    private JSONArray normalizeSubscribers(JSONObject data) {
        Object d=data.opt("data");
        if(d instanceof JSONArray)return (JSONArray)d;
        if(d instanceof JSONObject){JSONObject o=(JSONObject)d; if(o.opt("subscribers") instanceof JSONArray)return o.optJSONArray("subscribers");}
        if(data.opt("subscribers") instanceof JSONArray)return data.optJSONArray("subscribers");
        return new JSONArray();
    }

    private JSONArray normalizeArray(JSONObject data, String... keys) {
        for(String k:keys){Object v=data.opt(k); if(v instanceof JSONArray)return (JSONArray)v;}
        return new JSONArray();
    }
'''
new = '''    private JSONArray normalizeSubscribers(JSONObject data) {
        if (data == null) return new JSONArray();
        Object root = data.opt("_root"); if (root instanceof JSONArray) return (JSONArray) root;
        Object d=data.opt("data");
        if(d instanceof JSONArray)return (JSONArray)d;
        if(d instanceof JSONObject){JSONObject o=(JSONObject)d; if(o.opt("subscribers") instanceof JSONArray)return o.optJSONArray("subscribers"); if(o.opt("data") instanceof JSONArray)return o.optJSONArray("data");}
        if(data.opt("subscribers") instanceof JSONArray)return data.optJSONArray("subscribers");
        Object result=data.opt("result");
        if(result instanceof JSONArray)return (JSONArray)result;
        if(result instanceof JSONObject){JSONObject o=(JSONObject)result; if(o.opt("subscribers") instanceof JSONArray)return o.optJSONArray("subscribers");}
        return new JSONArray();
    }

    private JSONArray normalizeArray(JSONObject data, String... keys) {
        if (data == null) return new JSONArray();
        Object root=data.opt("_root"); if(root instanceof JSONArray)return (JSONArray)root;
        for(String k:keys){
            Object v=data.opt(k);
            if(v instanceof JSONArray)return (JSONArray)v;
            if(v instanceof JSONObject){
                JSONObject o=(JSONObject)v;
                for(String nested:keys){Object nv=o.opt(nested);if(nv instanceof JSONArray)return(JSONArray)nv;}
            }
        }
        Object d=data.opt("data"); if(d instanceof JSONArray)return(JSONArray)d;
        return new JSONArray();
    }
'''
if old not in s: raise SystemExit('v11 patch: normalize marker missing')
s = s.replace(old, new)

marker = '''    private String firstString(JSONObject o,String...keys){for(String k:keys){String v=o.optString(k,"");if(!v.isEmpty())return v;}return "";}'''
recent = '''    private JSONArray recentTokens() {
        SharedPreferences p = getSharedPreferences("gpt_manager", MODE_PRIVATE);
        String raw = p.getString("recent_tokens_v1", "[]");
        try { return new JSONArray(raw == null ? "[]" : raw); } catch (JSONException e) { return new JSONArray(); }
    }

    private void saveRecentToken(String token, String role, String label) {
        if (token == null || token.trim().length() < 4) return;
        token = token.trim();
        JSONArray old = recentTokens();
        JSONArray next = new JSONArray();
        try {
            next.put(new JSONObject().put("token", token).put("role", role).put("label", label == null ? "" : label).put("lastUsed", System.currentTimeMillis()));
            for (int i=0;i<old.length() && next.length()<10;i++) {
                JSONObject row=old.optJSONObject(i); if(row==null)continue;
                if(!token.equals(row.optString("token"))) next.put(row);
            }
            getSharedPreferences("gpt_manager", MODE_PRIVATE).edit().putString("recent_tokens_v1", next.toString()).apply();
        } catch (JSONException ignored) {}
    }

    private JSONObject recentTokensResult() {
        try { return ok().put("data", recentTokens()); } catch (JSONException e) { return error("خواندن توکن‌های اخیر ناموفق بود."); }
    }

    private JSONObject forgetRecentToken(String token) {
        JSONArray old=recentTokens(), next=new JSONArray();
        for(int i=0;i<old.length();i++){JSONObject row=old.optJSONObject(i);if(row!=null && !token.equals(row.optString("token")))next.put(row);}
        getSharedPreferences("gpt_manager", MODE_PRIVATE).edit().putString("recent_tokens_v1", next.toString()).apply();
        return ok();
    }

    private JSONObject clearRecentTokens() {
        getSharedPreferences("gpt_manager", MODE_PRIVATE).edit().remove("recent_tokens_v1").apply();
        return ok();
    }

''' + marker
if marker not in s: raise SystemExit('v11 patch: helper marker missing')
s = s.replace(marker, recent)

old = '''        setContentView(root);
        webView.loadUrl(DASHBOARD);'''
new = '''        setContentView(root);
        clearProxy();
        webView.loadUrl(DASHBOARD);'''
if old not in s: raise SystemExit('v11 patch: onCreate marker missing')
s = s.replace(old, new)

marker = '''    private String hostOf(String url) {'''
new = '''    @Override protected void onDestroy() {
        clearProxy();
        try { proxyExecutor.shutdownNow(); } catch (Throwable ignored) {}
        super.onDestroy();
    }

''' + marker
if marker not in s: raise SystemExit('v11 patch: host marker missing')
s = s.replace(marker, new)

required = [
    'case "recentTokens"',
    'trimmed.startsWith("[")',
    'recent_tokens_v1',
    'subscriptionLoadOk',
    'webView.removeJavascriptInterface("AndroidApi")',
]
for item in required:
    if item not in s: raise SystemExit('v11 patch missing: ' + item)
java.write_text(s, encoding='utf-8')

index.write_text((root / 'v11_index.html').read_text(encoding='utf-8'), encoding='utf-8')
icon.parent.mkdir(parents=True, exist_ok=True)
icon.write_text((root / 'v11_icon.xml').read_text(encoding='utf-8'), encoding='utf-8')

m = manifest.read_text(encoding='utf-8')
if 'android:icon="@drawable/ic_launcher"' not in m:
    m = m.replace('android:theme="@style/AppTheme"\n        android:label="GPT Manager"', 'android:theme="@style/AppTheme"\n        android:label="GPT Manager"\n        android:icon="@drawable/ic_launcher"\n        android:roundIcon="@drawable/ic_launcher"')
manifest.write_text(m, encoding='utf-8')

b = build.read_text(encoding='utf-8')
b = b.replace("versionCode 1\n        versionName '1.0.0'", "versionCode 2\n        versionName '1.1.0'")
build.write_text(b, encoding='utf-8')

print('Applied GPT Manager Android v1.1.0 patch')
