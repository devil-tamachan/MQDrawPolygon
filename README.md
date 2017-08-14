# MQDrawPolygon
ラスター画像をポリゴン化。メタセコイア４用プラグイン

## 更新履歴
 - v1.2 (20170814): 面に貼ったテクスチャ画像からの変換、立体テクスチャ機能追加
 - v1.1 (20170510): 面貼り追加

## ダウンロード
 - 64ビット版: https://github.com/devil-tamachan/MQDrawPolygon/releases
 - 32ビット版: 今のところありません。自分でビルドしてください

## インストール方法
 - MQDrawPolygon.dll, DrawPolygon.exe, CGAL-vc100-mt-4.9.1.dl_, libgmp-10.dl_, libmpfr-4.dl_ の５ファイルを　Metasequoia4_x64\Plugins\Select　フォルダへコピー、貼り付け
 - MQShrinkWrapを使っている場合、CGAL-vc100-mt-4.9.1.dl_ はスキップか上書きしてください
 - v1.1でlibgmp-10.dl_が追加されました。忘れずにコピー
 - v1.2でDrawPolygon.exe, libmpfr-4.dl_が追加されました。忘れずにコピー

## 使い方（面テクスチャ　高速/低速モード）
<img src="https://raw.githubusercontent.com/devil-tamachan/MQDrawPolygon/master/sample/rittai.png" />
<img src="https://raw.githubusercontent.com/devil-tamachan/MQDrawPolygon/master/sample/furiru2.png" />
<img src="https://raw.githubusercontent.com/devil-tamachan/MQDrawPolygon/master/sample/pakaan.png" />
 - 変換サンプル: sample/spikebelt.mqo, sample/furiru2.mqo
 - テクスチャに模様を書いて面に貼る
 - 変換したい面を選択
 - メタセコイア４を起動
 - メニュー、選択部処理、DrawPolygon 
 - 設定してOK (しばらく待つ)
 - 変換された新しいオブジェクトが追加
 - （高速モードの場合のみ、キャッシュオブジェクトがテクスチャ画像ごとに追加される）

## 使い方（クリップボードから入力）
 - SAI/AzPainter2/IrfanViewなどで画像を開く、または描く
 - 画像をコピー
 - メタセコイア４を起動
 - メニュー、選択部処理、DrawPolygon 
 - 設定してOK (しばらく待つ)
 - 変換された新しいオブジェクトが追加
 - 面を貼る （現バージョンは線のみ出力。面は貼られません）

## 注意事項
 - 設定によってはめちゃくちゃ時間がかかります。最悪メタセコイア４を強制終了できるように、新しくメタセコイア４を起動することを推奨します
 - ライン抽出: "無し (中央線)"を選ぶ場合、できるだけ細い線(1px推奨)で描いてください。線がつながっている必要はありません。太い線だとめちゃくちゃ時間がかかります。AzDrawing2のドットペン推奨
 - 出力頂点数はできるだけ少なくしてください。カクカクしていたら徐々に増やしてみてください

## 変換サンプル
<img src="https://github.com/devil-tamachan/MQDrawPolygon/raw/master/sample/ribbon_result.png" />
<img src="https://github.com/devil-tamachan/MQDrawPolygon/raw/master/sample/color1_result.png" />
<img src="https://github.com/devil-tamachan/MQDrawPolygon/raw/master/sample/centerline1_result.png" />

## メタセコイア４ プラグインリンク集
 - https://docs.google.com/spreadsheets/d/1HNXBzK-aAXoJ_Dp4ijxiaHkrvYuYB7sx4kyf-tdV0Hk/edit#gid=0
