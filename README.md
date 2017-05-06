# MQDrawPolygon
ラスター画像をポリゴン化。メタセコイア４用プラグイン

## ダウンロード
 - https://github.com/devil-tamachan/MQDrawPolygon/releases
 - 32ビット版は今のところありません。自分でビルドしてください

## インストール方法
 - MQDrawPolygon.dll, CGAL-vc100-mt-4.9.1.dl_ の２ファイルを　Metasequoia4_x64\Plugins\Select　フォルダへコピー、貼り付け
 - MQShrinkWrapを使っている場合、CGAL-vc100-mt-4.9.1.dl_ はスキップか上書きしてください

## 使い方
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
<img src="https://github.com/devil-tamachan/MQDrawPolygon/raw/master/sample/color1_result.png" />
<img src="https://github.com/devil-tamachan/MQDrawPolygon/raw/master/sample/centerline1_result.png" />
