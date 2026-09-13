from pathlib import Path
import base64
root=Path(__file__).resolve().parent
raw=base64.b64decode((root/'app_icon.b64').read_text().strip())
if not (raw.startswith(b'RIFF') and raw[8:12] == b'WEBP'):
    raise SystemExit('app icon is not WebP')
for rel in ['app/src/main/res/drawable/app_icon.webp','app/src/main/assets/app_icon.webp']:
    p=root/rel; p.parent.mkdir(parents=True,exist_ok=True); p.write_bytes(raw)
print(f'Prepared app icon: {len(raw)} bytes')
