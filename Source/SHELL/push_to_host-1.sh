#!/bin/bash
# schのソースや計算条件、計算結果などをhost-1へとコピーする
# OUTやSAVEは非常に容量が大きいこと、コピーに時間をとることからコピーはしない
# 引数に、その計算条件を表現するようなディレクトリ名を入れることを想定している
#  i.e. ./push_to_host-1.sh Case105_m1q1B1u1v1E1e-6kE1
# のようにする．Case番号さえ間違っていなければ問題ない．
#
# 具体的には以下のことをする．
# 1. コピーする前に，結果からgifアニメやlocus_info.epsを作っておく．
# 2. このホスト上に今回の計算結果やソースコードを新たにディレクトリを作って保存
# 3. OUT, SAVE以外をhost-1へコピー
#
# SAVEの容量がそれなりに大きいので(4Gくらい)，計算結果として格納するのを止める．

set -eu # -e: コマンド失敗があればスクリプトを打ち止め，-u: 未定義の変数があると打ち止め

# 使い方の説明
function Usage(){
	echo "Usage:"
	echo "	./push_to_host-1.sh [Directory name]"
	echo "	For example:"
	echo "	./push_to_host-1.sh Case105_m1q1B1u1v1E1e-6kE1"
	exit 2;
}

if [ $# -ne 1 ]; then
	Usage
	exit 2
fi

# Directory name and destination in host-1
# CASEは引数で指定した解析条件のCase番号と条件の羅列
CASE=$1
# 自分自身のパスを取得
PWD=`pwd` # /home/wataru_k/Sinusoidal/sch_code/V_cos/V_cos_ver1.0/SHELL とかのはず
# PWDからSHELLの文字列を取り除く
SRC="${PWD%%/SHELL}"
# 計算結果を格納するディレクトリを指定
#DATA="$HOME/Sinusoidal/schdata/u0_eq_0" # 自分のホームディレクトリを表す"~"を使うとrsyncでエラーが起きる
#DATA="$HOME/Sinusoidal/schdata/u0_neq_0" # 自分のホームディレクトリを表す"~"を使うとrsyncでエラーが起きる
DATA="$HOME/High_Accuracy_Sinusoidal/sch_data" # 自分のホームディレクトリを表す"~"を使うとrsyncでエラーが起きる
# 格納されるパスを生成
DST="${DATA}/${CASE}"
#echo "Debug: CASE = $CASE"
#echo "Debug: PWD = $PWD"
#echo "Debug: SRC = $SRC"
#echo "Debug: DATA = $DATA"
#echo "Debug: DST = $DST"
# 今のhostでもhost-1でもデータを格納するディレクトリは同じパスでかけるとする

## ディレクトリ名の確認
echo "********************************************************************************"
echo "	CASE = $CASE"
echo "	DST  = $DST"
echo "********************************************************************************"
echo "OK? [Y/n]"
read ANSWER
case $ANSWER in
	"" | "Y" | "y" | "yes" | "Yes" | "YES" )
	  	echo "OK. Processing forward!"
		;;
	*)
	  	echo "Process is canceled. So, any files is not moved."
		exit 1
		;;
esac

#echo "exit 0 for debug."; exit 0; # for debug

# renew locus_info.eps
echo "--- Drawing figures ---"
cd "${SRC}/PLT"
echo "gnuplot toEPS_multi3.plt"
gnuplot toEPS_multi3.plt

# gifファイルが存在していないか，
#GIF="$SRC/PLT/real8.gif"
#OUTlast=`ls ../OUT | tail -n1`
# 既にあるgifファイルのタイムスタンプがOUT内の最後の.rhoファイルより過去の場合
# if [ ! -e "$GIF" ] || [ "$GIF" -ot "$OUTlast" ] ; then
# 	echo "gnuplot bb8.plt ... This may take some time."
# 	gnuplot bb8.plt
# 	echo "gif animation done!"
# fi
# -> 時間がかかるので自動化は止める．

cd $SRC # sch.cuや000.datがあるディレクトリへ
echo "--- Copying data on this host ---"
echo "mkdir -p ${DST}"
mkdir -p ${DST}
# exit 0; # for debug

echo "cp and mv"
#cp Makefile *.cu *.c sch.inp potential.dat DEHINT.* ${DST}
cp Makefile *.cu *.c sch.inp DEHINT.* ${DST}
cp -r PLT ${DST}
cp -r SHELL ${DST}
mv aaa.out 000.dat ${DST}
# mv OUT SAVE  ${DST}
mv OUT ${DST}; echo "mv OUT"
mkdir OUT; echo "mkdir OUT"
# SAVEは保管しない
rm SAVE/*.bin; echo "rm SAVE/*.bin"

echo "--- Copying data to host-1 ---"
# どのホストで計算したかをコピーするファイルに含めておく
echo "touch \"${DST}/from_`hostname`\""
touch "${DST}/from_`hostname`"
#echo "aaa.out 000.dat ... etc are put into ${DST}, except OUT and SAVE"
echo "rysnc to host-1"
# シェルスクリプト内で転送元のファイルのパスを書くときには"~"を使うと思ったように展開されないので使用注意!
#echo "rsync -rtunv --exclude="OUT" --exclude="SAVE" ${DST}/ host-1.qe.eng.hokudai.ac.jp:${DST}"
#rsync -rtuvn --exclude="OUT" --exclude="SAVE" ${DST}/ host-1.qe.eng.hokudai.ac.jp:${DST}
echo "rsync -rtuv --exclude="OUT" ${DST}/ host-1.qe.eng.hokudai.ac.jp:${DST}"
      rsync -rtuv --exclude="OUT" ${DST}/ host-1.qe.eng.hokudai.ac.jp:${DST}
# rsyncのオプションの説明
# -r: 再帰的にコピー
# -u: 上書きしない
# -t: タイムスタンプそのまま
# -v: 冗長な出力．何が送られたとか詳しく表示する
# -n: 試験モード，実際にコピーはせずに何をするかを表示する
# 転送元のパスの最後に"/"がないとディレクトリをコピーするので，
# rsync aaa host-1:aaa
# とかすると，host-1上で，aaa/aaaというようにコピーされてしまう
# "/"をつけるとそのディレクトリ内のファイルらを転送先のディレクトリにコピーする．



