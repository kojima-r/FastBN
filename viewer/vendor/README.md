# viewer/vendor — 同梱している第三者ライブラリ

| ファイル | 内容 |
| --- | --- |
| `cosmos.gl-3.4.1.min.js` | [cosmos.gl](https://github.com/cosmosgl/graph) (`@cosmos.gl/graph@3.4.1`) の UMD ビルド (`dist/index.min.js` をそのままコピー)。依存 (luma.gl, d3, gl-matrix ほか) を含む単一ファイルで、`window.Cosmos.Graph` を公開する |
| `LICENSE.cosmos.gl.txt` | 上記のライセンス (MIT) |

`sha256(cosmos.gl-3.4.1.min.js) = d6343f78dda80667e92c6902918763e72a1b03b5ba165c1581f5b7adb7d06b60`

ビルド手順を持ち込まず**オフラインで動く**ようにするため同梱している (npm / bundler は不要)。
更新・再取得は次のいずれかで行う。

```bash
python3 viewer/serve.py --fetch-vendor              # 同梱版が無い / 壊れたときに取り直す
python3 viewer/serve.py --fetch-vendor --vendor-version 3.5.0   # 別のバージョンに上げる
```

手で取る場合:

```bash
curl -sSL https://registry.npmjs.org/@cosmos.gl/graph/-/graph-3.4.1.tgz | tar xz -O \
  package/dist/index.min.js > viewer/vendor/cosmos.gl-3.4.1.min.js
```

バージョンを変えたら `viewer/index.html` の `<script src=...>` も合わせて更新する
(`serve.py --fetch-vendor` は自動で書き換える)。
