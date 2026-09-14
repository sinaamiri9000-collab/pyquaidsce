from pathlib import Path

root = Path('.')
java = root / 'app/src/main/java/com/gptyar/manager/MainActivity.java'
index = root / 'app/src/main/assets/index.html'
build = root / 'app/build.gradle'

s = java.read_text(encoding='utf-8')
s = s.replace('GPTManagerAndroid/1.1', 'GPTManagerAndroid/1.2')

old = '''                    case "clearRecentTokens": return clearRecentTokens().toString();
                    default: return error("درخواست ناشناخته است.").toString();'''
new = '''                    case "clearRecentTokens": return clearRecentTokens().toString();
                    case "adminNotifications": return adminNotifications().toString();
                    case "adminNotificationCreate": return adminNotificationCreate(body).toString();
                    case "adminNotificationUpdate": return adminNotificationUpdate(body).toString();
                    case "adminNotificationDelete": return adminNotificationDelete(body.optLong("notificationId", -1)).toString();
                    default: return error("درخواست ناشناخته است.").toString();'''
if old not in s:
    raise SystemExit('v12 patch: bridge marker missing')
s = s.replace(old, new)

# Match the original extension: proxy-info is fetched directly with token only.
s = s.replace(
    'JSONObject pb = new JSONObject().put("token", token).put("subscriptionId", subscriptionId);',
    'JSONObject pb = new JSONObject().put("token", token);'
)

old = '''                ProxySpec spec = ProxySpec.from(p);
                if (spec.host.isEmpty() || spec.port <= 0) throw new Exception("Proxy نامعتبر است.");
                applyWebViewProxy(spec);
                synchronized (stateLock) { currentProxyUser = spec.user; currentProxyPass = spec.pass; }
                step(steps, "proxy", true, "Proxy #" + (i+1) + " → " + spec.host + ":" + spec.port);

                JSONObject pingBody = new JSONObject().put("token", token).put("device_id", deviceId())
                        .put("proxy_id", p.optLong("PROXY_ID", p.optLong("proxy_id", -1))).put("isChecked", true);
                ApiResponse ping = api("/clients/proxy/ping", "POST", token, pingBody, spec, true);
                if (!ping.ok) ping = api("/clients/proxy/ping", "POST", token, pingBody, null, true);
                if (!ping.ok) throw new Exception("Proxy Auth: " + messageOf(ping));
                step(steps, "proxyAuth", true, "احراز هویت Proxy تکمیل شد.");'''
new = '''                ProxySpec spec = ProxySpec.from(p);
                if (spec.host.isEmpty() || spec.port <= 0) throw new Exception("Proxy نامعتبر است.");
                long proxyId = p.optLong("PROXY_ID", p.optLong("proxy_id", -1));
                if (proxyId <= 0) throw new Exception("شناسه Proxy نامعتبر است.");

                // The extension authorizes the selected proxy by calling the API DIRECTLY.
                // Routing this request through the proxy causes HTTPS CONNECT 403 before auth is registered.
                JSONObject pingBody = new JSONObject().put("token", token).put("device_id", deviceId())
                        .put("proxyId", proxyId)
                        .put("client", new JSONObject().put("product", "PLUGIN").put("version", VERSION).put("platform", "CHROME"));
                ApiResponse ping = null;
                for (int authTry=0; authTry<3; authTry++) {
                    ping = api("/clients/proxy/ping", "POST", token, pingBody, null, true);
                    if (ping.ok) break;
                    if (authTry < 2) Thread.sleep(500L * (authTry + 1));
                }
                if (ping == null || !ping.ok) throw new Exception("Proxy Auth: " + (ping == null ? "no response" : messageOf(ping)));
                step(steps, "proxyAuth", true, "احراز هویت Proxy تکمیل شد.");

                applyWebViewProxy(spec);
                synchronized (stateLock) { currentProxyUser = spec.user; currentProxyPass = spec.pass; }
                step(steps, "proxy", true, "Proxy #" + (i+1) + " → " + spec.host + ":" + spec.port);'''
if old not in s:
    raise SystemExit('v12 patch: proxy auth block missing')
s = s.replace(old, new)

# Cookie/session API is also direct in the original extension; only target service traffic uses the selected proxy.
s = s.replace(
    'cookiesResponse = api("/cookie/fetch", "POST", token, cb, spec, true);',
    'cookiesResponse = api("/cookie/fetch", "POST", token, cb, null, true);'
)

old = '''            } catch (Throwable t) {
                last = t;
                step(steps, "attempt", false, "Proxy #"+(i+1)+": "+(t.getMessage()==null?t.toString():t.getMessage()));
            }'''
new = '''            } catch (Throwable t) {
                last = t;
                clearProxy();
                synchronized (stateLock) { currentProxyUser = ""; currentProxyPass = ""; }
                step(steps, "attempt", false, "Proxy #"+(i+1)+": "+(t.getMessage()==null?t.toString():t.getMessage()));
            }'''
if old not in s:
    raise SystemExit('v12 patch: activation catch marker missing')
s = s.replace(old, new)

s = s.replace(
    'JSONArray noteRows = notes.ok ? normalizeArray(notes.data, "notifications", "data", "result") : new JSONArray();',
    'JSONArray noteRows = notes.ok ? normalizeNotifications(notes.data) : new JSONArray();'
)

marker = '''    private JSONObject adminCreate(JSONObject body) throws Exception {
        String token; synchronized (stateLock) { token = sessionToken; }
        JSONObject b = new JSONObject().put("name", body.optString("name"));
        if (!body.optString("expiresAt").isEmpty()) b.put("expiresAt", body.optString("expiresAt"));
        ApiResponse r = api("/accupdator/subscribers/create", "POST", token, b, null, true);
        return r.ok ? ok().put("data", r.data) : upstreamError(r, "ایجاد کاربر ناموفق بود.");
    }
'''
addition = marker + '''
    private ApiResponse adminNotificationRequest(String method, String[] paths, JSONObject body) throws Exception {
        String token, role;
        synchronized (stateLock) { token = sessionToken; role = sessionRole; }
        if (!"admin".equals(role) || token.isEmpty()) {
            return new ApiResponse(false, 403, new JSONObject().put("message", "دسترسی مدیر لازم است."));
        }
        ApiResponse last = null;
        for (String path : paths) {
            try {
                ApiResponse r = api(path, method, token, body, null, true);
                last = r;
                if (r.ok) return r;
            } catch (Throwable t) {
                last = new ApiResponse(false, 0, new JSONObject().put("message", t.getMessage() == null ? t.toString() : t.getMessage()));
            }
        }
        return last == null ? new ApiResponse(false, 404, new JSONObject().put("message", "مسیر مدیریت اعلان روی سرور پیدا نشد.")) : last;
    }

    private JSONObject adminNotifications() throws Exception {
        ApiResponse r = adminNotificationRequest("GET", new String[]{
                "/accupdator/notifications/all",
                "/accupdator/notifications",
                "/accupdator/notifications/list",
                "/accupdator/notices/all",
                "/accupdator/notices"
        }, null);
        if (!r.ok) return upstreamError(r, "دریافت اعلان‌ها ناموفق بود.");
        return ok().put("data", normalizeNotifications(r.data));
    }

    private JSONObject adminNotificationCreate(JSONObject body) throws Exception {
        String title = body.optString("title", "").trim();
        String message = body.optString("message", "").trim();
        if (title.isEmpty() || message.isEmpty()) return error("عنوان و متن اعلان لازم است.");
        JSONObject b = new JSONObject().put("title", title).put("message", message).put("active", body.optBoolean("active", true));
        ApiResponse r = adminNotificationRequest("POST", new String[]{
                "/accupdator/notifications/create",
                "/accupdator/notifications",
                "/accupdator/notices/create",
                "/accupdator/notices"
        }, b);
        return r.ok ? ok().put("data", r.data) : upstreamError(r, "ثبت اعلان ناموفق بود.");
    }

    private JSONObject adminNotificationUpdate(JSONObject body) throws Exception {
        long id = body.optLong("notificationId", body.optLong("id", -1));
        if (id <= 0) return error("شناسه اعلان معتبر نیست.");
        String title = body.optString("title", "").trim();
        String message = body.optString("message", "").trim();
        JSONObject b = new JSONObject().put("notificationId", id).put("id", id)
                .put("title", title).put("message", message).put("active", body.optBoolean("active", true));
        ApiResponse r = adminNotificationRequest("PATCH", new String[]{
                "/accupdator/notifications/update",
                "/accupdator/notifications/" + id,
                "/accupdator/notices/update",
                "/accupdator/notices/" + id
        }, b);
        if (!r.ok && (r.status == 404 || r.status == 405)) {
            r = adminNotificationRequest("PUT", new String[]{
                    "/accupdator/notifications/" + id,
                    "/accupdator/notices/" + id
            }, b);
        }
        return r.ok ? ok().put("data", r.data) : upstreamError(r, "ویرایش اعلان ناموفق بود.");
    }

    private JSONObject adminNotificationDelete(long id) throws Exception {
        if (id <= 0) return error("شناسه اعلان معتبر نیست.");
        JSONObject b = new JSONObject().put("notificationId", id).put("id", id);
        ApiResponse r = adminNotificationRequest("DELETE", new String[]{
                "/accupdator/notifications/delete",
                "/accupdator/notifications/" + id,
                "/accupdator/notices/delete",
                "/accupdator/notices/" + id
        }, b);
        return r.ok ? ok() : upstreamError(r, "حذف اعلان ناموفق بود.");
    }
'''
if marker not in s:
    raise SystemExit('v12 patch: adminCreate marker missing')
s = s.replace(marker, addition)

norm_marker = '''    private JSONArray normalizeSubscribers(JSONObject data) {'''
norm = '''    private JSONArray normalizeNotifications(JSONObject data) {
        JSONArray found = findNotificationArray(data, 0);
        return found == null ? new JSONArray() : found;
    }

    private JSONArray findNotificationArray(Object node, int depth) {
        if (node == null || node == JSONObject.NULL || depth > 6) return null;
        if (node instanceof JSONArray) return (JSONArray) node;
        if (!(node instanceof JSONObject)) return null;
        JSONObject o = (JSONObject) node;
        String[] keys = {"notifications", "active", "notices", "alerts", "data", "result", "items"};
        for (String k : keys) {
            Object v = o.opt(k);
            if (v instanceof JSONArray) return (JSONArray) v;
        }
        for (String k : keys) {
            Object v = o.opt(k);
            if (v instanceof JSONObject) {
                JSONArray a = findNotificationArray(v, depth + 1);
                if (a != null) return a;
            }
        }
        if (o.has("title") || o.has("message") || o.has("text") || o.has("body") || o.has("content")) {
            JSONArray one = new JSONArray(); one.put(o); return one;
        }
        return null;
    }

''' + norm_marker
if norm_marker not in s:
    raise SystemExit('v12 patch: normalize marker missing')
s = s.replace(norm_marker, norm)

required = [
    '.put("proxyId", proxyId)',
    'api("/clients/proxy/ping", "POST", token, pingBody, null, true)',
    'api("/cookie/fetch", "POST", token, cb, null, true)',
    'case "adminNotifications"',
    'normalizeNotifications(notes.data)',
]
for item in required:
    if item not in s:
        raise SystemExit('v12 patch missing marker: ' + item)

# Ensure the broken proxy-auth path is gone.
if 'api("/clients/proxy/ping", "POST", token, pingBody, spec, true)' in s:
    raise SystemExit('v12 patch: proxy ping is still routed through selected proxy')
if 'api("/cookie/fetch", "POST", token, cb, spec, true)' in s:
    raise SystemExit('v12 patch: cookie fetch is still routed through selected proxy')

java.write_text(s, encoding='utf-8')
index.write_text((root / 'v12_index.html').read_text(encoding='utf-8'), encoding='utf-8')

b = build.read_text(encoding='utf-8')
b = b.replace("versionCode 2\n        versionName '1.1.0'", "versionCode 3\n        versionName '1.2.0'")
if "versionName '1.2.0'" not in b:
    raise SystemExit('v12 patch: version update failed')
build.write_text(b, encoding='utf-8')

print('Applied GPT Manager Android v1.2.0 patch')
