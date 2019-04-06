#!/bin/bash
set -eu
#*****************************************************************************************
# schのソースや計算条件、計算結果などをhost-1へとコピーする
# OUTやSAVEは非常に容量が大きいこと、コピーに時間をとることからコピーはしない
# 引数にその計算条件を表現するディレクトリ名を入れる
#  i.e. ./push_to_host-1.sh Case105_m1q1B1u1v1E1e-6kE1
# author Wataru kosaka   revised by shimizu
#*****************************************************************************************

#*****************************************************************************************
# 内容
# 1. このホスト上に今回の計算結果やソースコードを新たにディレクトリを作って保存
# 2. OUT, SAVE以外をhost-1へコピー
#
# SAVEは容量が大きい(4G程度)ので格納しない．
#*****************************************************************************************
#使い方の説明
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
#*****************************************************************************************
#デバッグプリント debug_print()
debug_print(){
	echo "End of this program"
	exit 0;
}
#*****************************************************************************************
# Directory name and destination in host-1

CASE=$1 										# 引数をcaseに格納
PWD=`pwd` 										# /home/wataru_k/Sinusoidal/sch_code/V_cos/V_cos_ver1.0/SHELL とかのはず
SOURCE="${PWD%%/SHELL}"							# PWDからSHELLの文字列を取り除く
DATA="$HOME/Landau_E=y^4/sch_data" 				# 自分のホームディレクトリを表す"~"を使うとrsyncでエラーが起きる
DESTINATION="${DATA}/${CASE}"							# 格納されるパスを生成
#---------------------------------------------------------------------------------------
# ディレクトリ名の確認
#---------------------------------------------------------------------------------------
echo "********************************************************************************"
echo "	CASE = $CASE"
echo "	destination  = $DESTINATION"
echo "********************************************************************************"
echo "OK? [Y/n]"
read ANSWER
case $ANSWER in
	"" | "Y" | "y" | "yes" | "Yes" | "YES" )
	  	echo "OK. Start copying now!"
		;;
	*)
	  	echo "The process is canceled. No files are moved."
		exit 1
		;;
esac
#***************************************************************************************
#	コピー処理開始
#***************************************************************************************
#---------------------------------------------------------------------------------------
#	ホスト内へバックアップをコピー
#---------------------------------------------------------------------------------------
cd $SOURCE 															# sch.cuや000.datがあるディレクトリへ
echo "--- Copying data on the home directory of this host ---"
mkdir -p ${DESTINATION}
cp Makefile *.cu *.c sch.inp DEHINT.* ${DESTINATION}
cp -r PLT ${DESTINATION}
cp -r SHELL ${DESTINATION}
mv aaa.out 000.dat ${DESTINATION}
# mv OUT SAVE  ${DESTINATION}
mv OUT ${DESTINATION}; echo "mv OUT"
mkdir OUT; echo "mkdir OUT"
rm SAVE/*.bin; echo "rm SAVE/*.bin"
#---------------------------------------------------------------------------------------
#	host-1へコピー
#---------------------------------------------------------------------------------------
echo "--- Copying data to host-1 ---"
touch "${DESTINATION}/from_`hostname`"
rsync -rtuv --exclude="OUT" ${DESTINATION}/ host-1.qe.eng.hokudai.ac.jp:${DESTINATION}