#!/bin/sh
#*************************************************************************
# execute.sh データファイルをエクセルに自動で取り込むためのプログラム群を
#				 execute（実行する）シェルスクリプト shimizu 2019/01/12
#*************************************************************************

#*************************************************************************
#内容
# 1.スクリプト群を実行するディレクトリ（一個上）へ一旦コピーする(4.で消す)
# 2.データファイル(tsv)をcsvを経由してxlsxファイルへ変換
# 3.xlsxファイルをテンプレートへ書き込む
# 4.いらないファイルを消去
#*************************************************************************

########################################################################1.
cp copy_sheet.py ../
cp copy_text_to_xlsx.py ../
cp CsvtoXlsx.py ../
cp TsvtoCsv.py ../
cd ../
echo 'Coping neccessary files has been completed.'
#文頭のコメント記号で要素化しないようにスペースを消去
sed -i "1 s/#    time/#____time/g" 000.dat
sed -i "1 s/#    time/#____time/g" 000.ufm

########################################################################2.
python3 TsvtoCsv.py 000.dat
python3 TsvtoCsv.py 000.ufm
python3 TsvtoCsv.py sch.inp
python3 CsvtoXlsx.py 000_dat.csv float
python3 CsvtoXlsx.py 000_ufm.csv float
python3 CsvtoXlsx.py sch_inp.csv
########################################################################3.
echo 'get started to copy 000.dat to Landau.xlsx ...'
python3 copy_sheet.py 000_dat.xlsx Sheet1 Landau_LE4_template_shimizu.xlsx 000.dat
echo 'get started to copy 000.ufm to Landau.xlsx ...'
python3 copy_sheet.py 000_ufm.xlsx Sheet1 Landau_LE4_template_shimizu.xlsx 000.ufm
echo 'get started to copy sch.inp to Landau.xlsx ...'
python3 copy_sheet.py sch_inp.xlsx Sheet1 Landau_LE4_template_shimizu.xlsx sch.inp
echo 'get started to copy potential.cu to Landau.xlsx ...'
python3 copy_text_to_xlsx.py potential.cu Landau_LE4_template_shimizu.xlsx potential.cu

echo 'Coping has been completed.'
########################################################################4.
rm copy_sheet.py copy_text_to_xlsx.py CsvtoXlsx.py TsvtoCsv.py
rm sch_inp.xlsx 000_dat.xlsx 000_ufm.xlsx 000_dat.csv 000_ufm.csv sch_inp.csv
echo 'unnecessaty files have been deleted.'
