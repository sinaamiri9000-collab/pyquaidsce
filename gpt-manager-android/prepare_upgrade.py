from pathlib import Path
import base64, gzip
root=Path(__file__).resolve().parent
packed=base64.b64decode((root/'upgrade_v11.patch.gz.b64').read_text().strip())
patch=gzip.decompress(packed)
(root/'upgrade_v11.patch').write_bytes(patch)
print(f'Prepared source patch: {len(patch)} bytes')
