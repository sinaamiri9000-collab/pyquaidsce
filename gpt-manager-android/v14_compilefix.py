from pathlib import Path
p=Path('app/src/main/java/com/gptyar/manager/MainActivity.java')
s=p.read_text(encoding='utf-8')
old='''    @Override protected void onDestroy() {
        clearProxy();
        try { proxyExecutor.shutdownNow(); } catch (Throwable ignored) {}
        super.onDestroy();
    }

'''
if old in s:
    s=s.replace(old,'',1)
# Assert exactly one lifecycle method remains.
if s.count('protected void onDestroy()') != 1:
    raise SystemExit(f'expected exactly one onDestroy, found {s.count("protected void onDestroy()")}')
p.write_text(s,encoding='utf-8')
print('Removed duplicate legacy onDestroy; exactly one remains')
