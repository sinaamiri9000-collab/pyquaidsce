package com.gptyar.manager;

import android.annotation.SuppressLint;
import android.app.Activity;
import android.content.ClipData;
import android.content.ClipboardManager;
import android.content.Context;
import android.content.SharedPreferences;
import android.graphics.Color;
import android.net.http.SslError;
import android.os.Bundle;
import android.os.Handler;
import android.os.Looper;
import android.view.Gravity;
import android.view.View;
import android.webkit.CookieManager;
import android.webkit.HttpAuthHandler;
import android.webkit.JavascriptInterface;
import android.webkit.SslErrorHandler;
import android.webkit.WebResourceRequest;
import android.webkit.WebSettings;
import android.webkit.WebView;
import android.webkit.WebViewClient;
import android.widget.Button;
import android.widget.LinearLayout;
import android.widget.TextView;
import android.widget.Toast;

import androidx.webkit.ProxyConfig;
import androidx.webkit.ProxyController;
import androidx.webkit.WebViewFeature;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.InetSocketAddress;
import java.net.Proxy;
import java.net.URL;
import java.nio.charset.StandardCharsets;
import java.util.Base64;
import java.util.Locale;
import java.util.UUID;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

public class MainActivity extends Activity {
    private static final String API = "https://api.gptyar.com";
    private static final String VERSION = "4.6.3";
    private static final String DASHBOARD = "file:///android_asset/index.html";

    private WebView webView;
    private LinearLayout toolbar;
    private TextView toolbarTitle;
    private final Handler main = new Handler(Looper.getMainLooper());
    private final ExecutorService proxyExecutor = Executors.newSingleThreadExecutor();
    private final Object stateLock = new Object();

    private String sessionToken = "";
    private String sessionRole = "";
    private String currentProxyUser = "";
    private String currentProxyPass = "";
    private JSONArray cachedSubscribers = new JSONArray();

    @Override
    @SuppressLint({"SetJavaScriptEnabled", "JavascriptInterface"})
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);

        LinearLayout root = new LinearLayout(this);
        root.setOrientation(LinearLayout.VERTICAL);
        root.setBackgroundColor(Color.rgb(7, 20, 38));

        toolbar = new LinearLayout(this);
        toolbar.setOrientation(LinearLayout.HORIZONTAL);
        toolbar.setGravity(Gravity.CENTER_VERTICAL);
        toolbar.setPadding(16, 8, 16, 8);
        toolbar.setBackgroundColor(Color.rgb(9, 28, 50));
        toolbar.setVisibility(View.GONE);

        Button back = new Button(this);
        back.setText("بازگشت به پنل");
        back.setOnClickListener(v -> returnToDashboard());
        toolbar.addView(back, new LinearLayout.LayoutParams(0, LinearLayout.LayoutParams.WRAP_CONTENT, 1));

        toolbarTitle = new TextView(this);
        toolbarTitle.setTextColor(Color.WHITE);
        toolbarTitle.setTextSize(16);
        toolbarTitle.setGravity(Gravity.END);
        toolbar.addView(toolbarTitle, new LinearLayout.LayoutParams(0, LinearLayout.LayoutParams.MATCH_PARENT, 2));

        webView = new WebView(this);
        WebSettings s = webView.getSettings();
        s.setJavaScriptEnabled(true);
        s.setDomStorageEnabled(true);
        s.setDatabaseEnabled(true);
        s.setAllowFileAccess(true);
        s.setAllowContentAccess(false);
        s.setMixedContentMode(WebSettings.MIXED_CONTENT_NEVER_ALLOW);
        s.setUserAgentString(s.getUserAgentString() + " GPTManagerAndroid/1.0");
        CookieManager.getInstance().setAcceptCookie(true);
        CookieManager.getInstance().setAcceptThirdPartyCookies(webView, true);
        webView.addJavascriptInterface(new NativeBridge(), "AndroidApi");
        webView.setWebViewClient(new WebViewClient() {
            @Override public boolean shouldOverrideUrlLoading(WebView view, WebResourceRequest request) { return false; }
            @Override public void onPageFinished(WebView view, String url) {
                boolean local = url != null && url.startsWith("file:///android_asset/");
                toolbar.setVisibility(local ? View.GONE : View.VISIBLE);
                toolbarTitle.setText(local ? "" : hostOf(url));
            }
            @Override public void onReceivedHttpAuthRequest(WebView view, HttpAuthHandler handler, String host, String realm) {
                String u, p;
                synchronized (stateLock) { u = currentProxyUser; p = currentProxyPass; }
                if (!u.isEmpty()) handler.proceed(u, p); else handler.cancel();
            }
            @Override public void onReceivedSslError(WebView view, SslErrorHandler handler, SslError error) {
                handler.cancel();
            }
        });

        root.addView(toolbar, new LinearLayout.LayoutParams(LinearLayout.LayoutParams.MATCH_PARENT, LinearLayout.LayoutParams.WRAP_CONTENT));
        root.addView(webView, new LinearLayout.LayoutParams(LinearLayout.LayoutParams.MATCH_PARENT, 0, 1));
        setContentView(root);
        webView.loadUrl(DASHBOARD);
    }

    private String hostOf(String url) {
        try { return new URL(url).getHost(); } catch (Exception e) { return "سرویس"; }
    }

    private void returnToDashboard() { webView.loadUrl(DASHBOARD); }

    @Override public void onBackPressed() {
        String u = webView.getUrl();
        if (u != null && !u.startsWith("file:///android_asset/")) { returnToDashboard(); return; }
        if (webView.canGoBack()) webView.goBack(); else super.onBackPressed();
    }

    private String deviceId() {
        SharedPreferences p = getSharedPreferences("gpt_manager", MODE_PRIVATE);
        String id = p.getString("device_id", "");
        if (id == null || id.length() < 8) {
            id = UUID.randomUUID().toString();
            p.edit().putString("device_id", id).apply();
        }
        return id;
    }

    private final class NativeBridge {
        @JavascriptInterface public String request(String action, String payload) {
            try {
                JSONObject body = payload == null || payload.isEmpty() ? new JSONObject() : new JSONObject(payload);
                switch (action) {
                    case "login": return login(body.optString("token")).toString();
                    case "logout": return logout().toString();
                    case "overview": return userOverview().toString();
                    case "activate": return activateCurrentUser(body.optLong("subscriptionId", -1)).toString();
                    case "adminSubscribers": return adminSubscribers().toString();
                    case "subscriberToken": return subscriberToken(body.optLong("subscriberId", -1)).toString();
                    case "adminActivate": return adminActivate(body.optLong("subscriberId", -1), body.optLong("targetSubscriptionId", -1)).toString();
                    case "adminUpdate": return adminUpdate(body).toString();
                    case "adminDelete": return adminDelete(body.optLong("subscriptionId", -1)).toString();
                    case "adminCreate": return adminCreate(body).toString();
                    default: return error("درخواست ناشناخته است.").toString();
                }
            } catch (Throwable t) { return error(t.getMessage() == null ? t.toString() : t.getMessage()).toString(); }
        }
        @JavascriptInterface public String getDeviceId() { return deviceId(); }
        @JavascriptInterface public void copy(String text) {
            main.post(() -> {
                ClipboardManager cb = (ClipboardManager) getSystemService(CLIPBOARD_SERVICE);
                cb.setPrimaryClip(ClipData.newPlainText("token", text == null ? "" : text));
                Toast.makeText(MainActivity.this, "کپی شد", Toast.LENGTH_SHORT).show();
            });
        }
        @JavascriptInterface public void openDashboard() { main.post(MainActivity.this::returnToDashboard); }
    }

    private JSONObject login(String token) throws Exception {
        token = token == null ? "" : token.trim();
        if (token.length() < 4) return error("توکن معتبر وارد کنید.");
        if (token.startsWith("@")) {
            JSONObject b = new JSONObject().put("token", token).put("version", VERSION);
            ApiResponse r = api("/accupdator/validate", "POST", token, b, null, true);
            if (!r.ok || r.data.optString("message").isEmpty() || !r.data.has("account") || !r.data.has("agent")) return upstreamError(r, "توکن مدیریتی معتبر نیست.");
            synchronized (stateLock) { sessionToken = token; sessionRole = "admin"; }
            JSONObject subs = adminSubscribers();
            return ok().put("role", "admin").put("profile", r.data).put("subscribers", subs.opt("data")).put("deviceId", deviceId());
        }

        JSONObject b = new JSONObject()
                .put("token", token).put("device_id", deviceId()).put("version", VERSION)
                .put("client", new JSONObject().put("product", "PLUGIN").put("version", VERSION).put("platform", "CHROME"));
        ApiResponse v = api("/client/validate", "POST", token, b, null, true);
        JSONObject notice = v.data.optJSONObject("notice");
        if (notice != null) return error(notice.optString("message", notice.optString("text", "ورود توسط سرور متوقف شد."))).put("notice", notice);
        if (!v.ok || !v.data.optBoolean("valid", false)) return upstreamError(v, "توکن معتبر نیست.");
        synchronized (stateLock) { sessionToken = token; sessionRole = "user"; }
        JSONObject overview = userOverview();
        overview.put("role", "user").put("deviceId", deviceId());
        return overview;
    }

    private JSONObject logout() throws Exception {
        clearServiceState();
        synchronized (stateLock) { sessionToken = ""; sessionRole = ""; cachedSubscribers = new JSONArray(); }
        return ok();
    }

    private JSONObject userOverview() throws Exception {
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

    private JSONObject activateCurrentUser(long subscriptionId) throws Exception {
        String token;
        synchronized (stateLock) { token = sessionToken; }
        if (token.isEmpty() || subscriptionId <= 0) return error("اشتراک معتبر نیست.");
        return activateToken(token, subscriptionId);
    }

    private JSONObject adminSubscribers() throws Exception {
        String token, role;
        synchronized (stateLock) { token = sessionToken; role = sessionRole; }
        if (!"admin".equals(role)) return error("دسترسی مدیر لازم است.");
        ApiResponse r = api("/accupdator/subscribers/all", "GET", token, null, null, true);
        if (!r.ok) return upstreamError(r, "دریافت کاربران ناموفق بود.");
        JSONArray rows = normalizeSubscribers(r.data);
        synchronized (stateLock) { cachedSubscribers = rows; }
        JSONArray safe = new JSONArray();
        for (int i=0;i<rows.length();i++) {
            JSONObject row = rows.optJSONObject(i); if (row == null) continue;
            JSONObject p = new JSONObject();
            long id = subscriberId(row);
            p.put("subscriptionId", id)
             .put("name", firstString(row, "name", "username", "title"))
             .put("expiresAt", firstString(row, "expiresAt", "expirationDate", "expires_at", "expiryDate"))
             .put("agentType", firstString(row, "agentType", "type"))
             .put("tokenAvailable", !findToken(row).isEmpty());
            safe.put(p);
        }
        return ok().put("data", safe);
    }

    private JSONObject subscriberToken(long subscriberId) throws Exception {
        JSONObject row = findSubscriber(subscriberId, true);
        if (row == null) return error("کاربر پیدا نشد.");
        String token = findToken(row);
        if (token.isEmpty()) return error("API برای این کاربر توکن برنگرداند.");
        return ok().put("token", token);
    }

    private JSONObject adminActivate(long subscriberId, long targetSubscriptionId) throws Exception {
        JSONObject row = findSubscriber(subscriberId, true);
        if (row == null) return error("کاربر پیدا نشد.");
        String userToken = findToken(row);
        if (userToken.isEmpty()) return error("توکن کاربر در پاسخ API موجود نیست.");

        JSONObject vb = new JSONObject().put("token", userToken).put("device_id", deviceId()).put("version", VERSION)
                .put("client", new JSONObject().put("product", "PLUGIN").put("version", VERSION).put("platform", "CHROME"));
        ApiResponse vr = api("/client/validate", "POST", userToken, vb, null, true);
        if (!vr.ok || !vr.data.optBoolean("valid", false)) return upstreamError(vr, "توکن کاربر برای این دستگاه معتبر نشد.");

        ApiResponse ar = api("/clients/subscription/agents-with-subscriptions", "GET", userToken, null, null, true);
        if (!ar.ok) return upstreamError(ar, "دریافت اشتراک‌های کاربر ناموفق بود.");
        JSONArray candidates = activationCandidates(normalizeArray(ar.data, "agents", "data"));
        if (targetSubscriptionId <= 0 && candidates.length() != 1) return ok().put("requiresSelection", true).put("candidates", candidates);
        long sid = targetSubscriptionId > 0 ? targetSubscriptionId : candidates.optJSONObject(0).optLong("subscriptionId", -1);
        return activateToken(userToken, sid);
    }

    private JSONObject activateToken(String token, long subscriptionId) throws Exception {
        JSONObject diag = new JSONObject().put("subscriptionId", subscriptionId).put("deviceId", deviceId()).put("steps", new JSONArray());
        JSONArray steps = diag.getJSONArray("steps");
        step(steps, "identity", true, "Device ID: " + deviceId());

        ApiResponse ar = api("/clients/subscription/agents-with-subscriptions", "GET", token, null, null, true);
        if (!ar.ok) return activationError(diag, "اشتراک‌ها دریافت نشدند: " + messageOf(ar));
        SubscriptionRef ref = findSubscription(normalizeArray(ar.data, "agents", "data"), subscriptionId);
        if (ref == null) return activationError(diag, "اشتراک برای این توکن پیدا نشد.");
        step(steps, "subscription", true, ref.type + " / " + ref.accountName);

        JSONObject pb = new JSONObject().put("token", token).put("subscriptionId", subscriptionId);
        ApiResponse pr = api("/client/proxy-info", "POST", "", pb, null, false);
        if (!pr.ok) return activationError(diag, "دریافت Proxy ناموفق بود: " + messageOf(pr));
        JSONArray proxies = normalizeArray(pr.data, "proxies", "data");
        if (proxies.length() == 0) return activationError(diag, "برای این اشتراک Proxy برگردانده نشد.");

        Throwable last = null;
        for (int i=0;i<Math.min(3, proxies.length());i++) {
            JSONObject p = proxies.optJSONObject(i); if (p == null) continue;
            try {
                ProxySpec spec = ProxySpec.from(p);
                if (spec.host.isEmpty() || spec.port <= 0) throw new Exception("Proxy نامعتبر است.");
                applyWebViewProxy(spec);
                synchronized (stateLock) { currentProxyUser = spec.user; currentProxyPass = spec.pass; }
                step(steps, "proxy", true, "Proxy #" + (i+1) + " → " + spec.host + ":" + spec.port);

                JSONObject pingBody = new JSONObject().put("token", token).put("device_id", deviceId())
                        .put("proxy_id", p.optLong("PROXY_ID", p.optLong("proxy_id", -1))).put("isChecked", true);
                ApiResponse ping = api("/clients/proxy/ping", "POST", token, pingBody, spec, true);
                if (!ping.ok) ping = api("/clients/proxy/ping", "POST", token, pingBody, null, true);
                if (!ping.ok) throw new Exception("Proxy Auth: " + messageOf(ping));
                step(steps, "proxyAuth", true, "احراز هویت Proxy تکمیل شد.");

                JSONObject cb = new JSONObject().put("token", token).put("device_id", deviceId()).put("subscriptionId", subscriptionId);
                ApiResponse cookiesResponse = null;
                JSONArray cookies = new JSONArray();
                for (int n=0;n<3;n++) {
                    cookiesResponse = api("/cookie/fetch", "POST", token, cb, spec, true);
                    cookies = collectCookies(cookiesResponse.data);
                    if (cookiesResponse.ok && cookies.length() > 0) break;
                    Thread.sleep(700L * (n+1));
                }
                if (cookiesResponse == null || !cookiesResponse.ok) throw new Exception("Session: " + (cookiesResponse == null ? "no response" : messageOf(cookiesResponse)));
                if (cookies.length() == 0) throw new Exception("No cookies returned");
                step(steps, "session", true, cookies.length() + " Cookie دریافت شد.");

                installCookies(cookies);
                step(steps, "cookies", true, "Cookieهای نشست اعمال شدند.");
                String target = serviceUrl(ref.type);
                main.post(() -> webView.loadUrl(target));
                return ok().put("activated", true).put("agentType", ref.type).put("subscriptionId", subscriptionId).put("diagnostic", diag);
            } catch (Throwable t) {
                last = t;
                step(steps, "attempt", false, "Proxy #"+(i+1)+": "+(t.getMessage()==null?t.toString():t.getMessage()));
            }
        }
        clearServiceState();
        return activationError(diag, last == null ? "فعال‌سازی ناموفق بود." : (last.getMessage()==null?last.toString():last.getMessage()));
    }

    private JSONObject adminUpdate(JSONObject body) throws Exception {
        String token; synchronized (stateLock) { token = sessionToken; }
        JSONObject b = new JSONObject().put("subscriptionId", body.optLong("subscriptionId")).put("expiresAt", body.optString("expiresAt"));
        ApiResponse r = api("/accupdator/subscribers/update", "PATCH", token, b, null, true);
        return r.ok ? ok().put("data", r.data) : upstreamError(r, "به‌روزرسانی ناموفق بود.");
    }
    private JSONObject adminDelete(long id) throws Exception {
        String token; synchronized (stateLock) { token = sessionToken; }
        ApiResponse r = api("/accupdator/subscribers/delete", "DELETE", token, new JSONObject().put("subscriptionId", id), null, true);
        return r.ok ? ok() : upstreamError(r, "حذف ناموفق بود.");
    }
    private JSONObject adminCreate(JSONObject body) throws Exception {
        String token; synchronized (stateLock) { token = sessionToken; }
        JSONObject b = new JSONObject().put("name", body.optString("name"));
        if (!body.optString("expiresAt").isEmpty()) b.put("expiresAt", body.optString("expiresAt"));
        ApiResponse r = api("/accupdator/subscribers/create", "POST", token, b, null, true);
        return r.ok ? ok().put("data", r.data) : upstreamError(r, "ایجاد کاربر ناموفق بود.");
    }

    private void applyWebViewProxy(ProxySpec p) throws Exception {
        if (!WebViewFeature.isFeatureSupported(WebViewFeature.PROXY_OVERRIDE)) throw new Exception("WebView این دستگاه Proxy Override را پشتیبانی نمی‌کند.");
        CountDownLatch latch = new CountDownLatch(1);
        ProxyConfig config = new ProxyConfig.Builder().addProxyRule(p.host + ":" + p.port).build();
        ProxyController.getInstance().setProxyOverride(config, proxyExecutor, latch::countDown);
        if (!latch.await(8, TimeUnit.SECONDS)) throw new Exception("زمان تنظیم Proxy تمام شد.");
    }

    private void clearProxy() {
        try {
            if (WebViewFeature.isFeatureSupported(WebViewFeature.PROXY_OVERRIDE)) {
                CountDownLatch latch = new CountDownLatch(1);
                ProxyController.getInstance().clearProxyOverride(proxyExecutor, latch::countDown);
                latch.await(4, TimeUnit.SECONDS);
            }
        } catch (Throwable ignored) {}
    }

    private void clearServiceState() {
        clearProxy();
        synchronized (stateLock) { currentProxyUser = ""; currentProxyPass = ""; }
        CountDownLatch latch = new CountDownLatch(1);
        main.post(() -> CookieManager.getInstance().removeAllCookies(v -> { CookieManager.getInstance().flush(); latch.countDown(); }));
        try { latch.await(3, TimeUnit.SECONDS); } catch (InterruptedException ignored) {}
    }

    private void installCookies(JSONArray cookies) throws Exception {
        CountDownLatch cleared = new CountDownLatch(1);
        main.post(() -> CookieManager.getInstance().removeAllCookies(v -> cleared.countDown()));
        cleared.await(5, TimeUnit.SECONDS);
        CountDownLatch done = new CountDownLatch(cookies.length());
        for (int i=0;i<cookies.length();i++) {
            JSONObject c = cookies.optJSONObject(i); if (c == null) { done.countDown(); continue; }
            String name = c.optString("name"), value = c.optString("value"), domain = c.optString("domain");
            if (name.isEmpty() || domain.isEmpty()) { done.countDown(); continue; }
            String url = c.optString("url");
            if (url.isEmpty()) url = "https://" + domain.replaceFirst("^\\.", "") + "/";
            StringBuilder line = new StringBuilder(name).append("=").append(value);
            line.append("; Path=").append(c.optString("path", "/"));
            if (!domain.isEmpty()) line.append("; Domain=").append(domain);
            if (c.optBoolean("secure", true)) line.append("; Secure");
            if (c.optBoolean("httpOnly", false)) line.append("; HttpOnly");
            String sameSite = c.optString("sameSite");
            if (!sameSite.isEmpty() && !"unspecified".equalsIgnoreCase(sameSite)) line.append("; SameSite=").append(capitalize(sameSite));
            final String u = url, cookieLine = line.toString();
            main.post(() -> CookieManager.getInstance().setCookie(u, cookieLine, accepted -> done.countDown()));
        }
        if (!done.await(12, TimeUnit.SECONDS)) throw new Exception("اعمال Cookieها کامل نشد.");
        main.post(() -> CookieManager.getInstance().flush());
    }

    private ApiResponse api(String path, String method, String token, JSONObject body, ProxySpec proxySpec, boolean auth) throws Exception {
        Proxy proxy = Proxy.NO_PROXY;
        if (proxySpec != null && !proxySpec.host.isEmpty() && proxySpec.port > 0) proxy = new Proxy(Proxy.Type.HTTP, new InetSocketAddress(proxySpec.host, proxySpec.port));
        HttpURLConnection c = (HttpURLConnection) new URL(API + path).openConnection(proxy);
        c.setRequestMethod(method);
        c.setConnectTimeout(12000); c.setReadTimeout(30000);
        c.setRequestProperty("Accept", "application/json");
        c.setRequestProperty("User-Agent", "GPTManagerAndroid/1.0 GPTYar/"+VERSION);
        if (auth && token != null && !token.isEmpty()) c.setRequestProperty("Authorization", "Bearer " + token);
        if (proxySpec != null && !proxySpec.user.isEmpty()) {
            String raw = proxySpec.user + ":" + proxySpec.pass;
            c.setRequestProperty("Proxy-Authorization", "Basic " + Base64.getEncoder().encodeToString(raw.getBytes(StandardCharsets.UTF_8)));
        }
        if (body != null) {
            byte[] bytes = body.toString().getBytes(StandardCharsets.UTF_8);
            c.setDoOutput(true); c.setRequestProperty("Content-Type", "application/json");
            try (OutputStream os = c.getOutputStream()) { os.write(bytes); }
        }
        int status = c.getResponseCode();
        InputStream is = status >= 200 && status < 400 ? c.getInputStream() : c.getErrorStream();
        String text = readAll(is);
        JSONObject data;
        try { data = text.isEmpty() ? new JSONObject() : new JSONObject(text); }
        catch (JSONException ex) { data = new JSONObject().put("message", text); }
        c.disconnect();
        return new ApiResponse(status >= 200 && status < 300, status, data);
    }

    private String readAll(InputStream is) throws Exception {
        if (is == null) return "";
        BufferedReader br = new BufferedReader(new InputStreamReader(is, StandardCharsets.UTF_8));
        StringBuilder sb = new StringBuilder(); String line;
        while ((line = br.readLine()) != null) sb.append(line);
        return sb.toString();
    }

    private JSONArray collectCookies(Object node) {
        JSONArray out = new JSONArray();
        collectCookiesRec(node, out, 0);
        return out;
    }
    private void collectCookiesRec(Object node, JSONArray out, int depth) {
        if (node == null || node == JSONObject.NULL || depth > 8) return;
        if (node instanceof JSONArray) {
            JSONArray a=(JSONArray)node; for(int i=0;i<a.length();i++) collectCookiesRec(a.opt(i),out,depth+1); return;
        }
        if (node instanceof JSONObject) {
            JSONObject o=(JSONObject)node;
            if (!o.optString("name").isEmpty() && o.has("value") && (!o.optString("domain").isEmpty() || !o.optString("url").isEmpty())) { out.put(o); return; }
            JSONArray names=o.names(); if(names!=null) for(int i=0;i<names.length();i++) collectCookiesRec(o.opt(names.optString(i)),out,depth+1);
        }
    }

    private JSONObject findSubscriber(long id, boolean refresh) throws Exception {
        if (refresh) adminSubscribers();
        synchronized (stateLock) {
            for(int i=0;i<cachedSubscribers.length();i++) {
                JSONObject r=cachedSubscribers.optJSONObject(i); if(r!=null && subscriberId(r)==id) return r;
            }
        }
        return null;
    }

    private long subscriberId(JSONObject r) {
        long v = r.optLong("subscriptionId", -1); if(v>0)return v;
        v=r.optLong("id",-1); if(v>0)return v;
        JSONObject s=r.optJSONObject("subscription"); return s==null?-1:s.optLong("id",-1);
    }

    private String findToken(Object node) { return findTokenRec(node,0); }
    private String findTokenRec(Object node,int depth) {
        if(!(node instanceof JSONObject)||depth>5)return "";
        JSONObject o=(JSONObject)node;
        String[] keys={"token","pluginToken","subscriberToken","userToken","clientToken","accessToken"};
        for(String k:keys){String v=o.optString(k,"").trim();if(v.length()>=6)return v;}
        String[] nested={"user","subscriber","subscription","account","client"};
        for(String k:nested){String v=findTokenRec(o.optJSONObject(k),depth+1);if(!v.isEmpty())return v;}
        return "";
    }

    private JSONArray normalizeSubscribers(JSONObject data) {
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

    private JSONArray activationCandidates(JSONArray agents) {
        JSONArray out=new JSONArray();
        for(int i=0;i<agents.length();i++){
            JSONObject a=agents.optJSONObject(i);if(a==null)continue;String type=agentType(a);
            JSONArray subs=subscriptions(a);
            for(int j=0;j<subs.length();j++){
                JSONObject s=subs.optJSONObject(j);if(s==null)continue;long id=s.optLong("id",s.optLong("subscriptionId",-1));
                if(id<=0)continue;
                JSONObject account=s.optJSONObject("account");
                try {
                    out.put(new JSONObject().put("subscriptionId",id).put("agentType",type)
                            .put("accountName",account==null?"":account.optString("name"))
                            .put("expiresAt",s.optString("expiresAt",s.optString("expirationDate",""))));
                } catch (JSONException ignored) {}
            }
        }
        return out;
    }

    private SubscriptionRef findSubscription(JSONArray agents,long id){
        for(int i=0;i<agents.length();i++){
            JSONObject a=agents.optJSONObject(i);if(a==null)continue;JSONArray subs=subscriptions(a);
            for(int j=0;j<subs.length();j++){
                JSONObject s=subs.optJSONObject(j);if(s==null)continue;long sid=s.optLong("id",s.optLong("subscriptionId",-1));
                if(sid==id){JSONObject ac=s.optJSONObject("account");return new SubscriptionRef(agentType(a), ac==null?"":ac.optString("name"),s,a);}
            }
        }return null;
    }
    private JSONArray subscriptions(JSONObject a){
        Object x=a.opt("userSubscriptions");if(x instanceof JSONArray)return(JSONArray)x;
        x=a.opt("subscriptions");if(x instanceof JSONArray)return(JSONArray)x;
        x=a.opt("user_subscriptions");if(x instanceof JSONArray)return(JSONArray)x;
        return new JSONArray();
    }
    private String agentType(JSONObject a){String x=firstString(a,"type","agentType","key","name");return x.toUpperCase(Locale.ROOT);}

    private String serviceUrl(String type){
        switch(type.toUpperCase(Locale.ROOT)){
            case "CLAUDE":return "https://claude.ai/";
            case "GROK":return "https://grok.com/";
            case "MIDJOURNEY":return "https://www.midjourney.com/";
            default:return "https://chatgpt.com/";
        }
    }

    private String firstString(JSONObject o,String...keys){for(String k:keys){String v=o.optString(k,"");if(!v.isEmpty())return v;}return "";}
    private String messageOf(ApiResponse r){return r.data.optString("message",r.data.optString("error","HTTP "+r.status));}
    private JSONObject upstreamError(ApiResponse r,String fallback){
        JSONObject o = error(r.data.optString("message",r.data.optString("error",fallback)));
        try { o.put("status",r.status).put("upstream",r.data); } catch (JSONException ignored) {}
        return o;
    }
    private JSONObject ok(){try{return new JSONObject().put("ok",true);}catch(Exception e){return new JSONObject();}}
    private JSONObject error(String m){try{return new JSONObject().put("ok",false).put("error",m==null?"خطای ناشناخته":m);}catch(Exception e){return new JSONObject();}}
    private JSONObject activationError(JSONObject diag,String m){try{diag.put("error",m);return error(m).put("diagnostic",diag);}catch(Exception e){return error(m);}}
    private void step(JSONArray a,String stage,boolean ok,String detail){
        try { a.put(new JSONObject().put("stage",stage).put("ok",ok).put("detail",detail)); } catch (JSONException ignored) {}
    }
    private String capitalize(String s){if(s==null||s.isEmpty())return s;return Character.toUpperCase(s.charAt(0))+s.substring(1).toLowerCase(Locale.ROOT);}

    private static final class ApiResponse { final boolean ok; final int status; final JSONObject data; ApiResponse(boolean o,int s,JSONObject d){ok=o;status=s;data=d;} }
    private static final class SubscriptionRef { final String type,accountName; final JSONObject subscription,agent; SubscriptionRef(String t,String n,JSONObject s,JSONObject a){type=t;accountName=n;subscription=s;agent=a;} }
    private static final class ProxySpec {
        final String host,user,pass; final int port;
        ProxySpec(String h,int p,String u,String pw){host=h;port=p;user=u;pass=pw;}
        static ProxySpec from(JSONObject o){
            String host=o.optString("PROXY_HOST",o.optString("host",o.optString("proxyHost","")));
            int port=o.optInt("PROXY_PORT",o.optInt("port",o.optInt("proxyPort",-1)));
            String server=o.optString("PROXY_SERVER",o.optString("server",""));
            if((host.isEmpty()||port<=0)&&server.contains(":")){
                int idx=server.lastIndexOf(':');host=server.substring(0,idx).replaceFirst("^https?://","");
                try{port=Integer.parseInt(server.substring(idx+1));}catch(Exception ignored){}
            }
            String user=o.optString("PROXY_USERNAME",o.optString("username",o.optString("user","")));
            String pass=o.optString("PROXY_PASSWORD",o.optString("password",o.optString("pass","")));
            return new ProxySpec(host,port,user,pass);
        }
    }
}
