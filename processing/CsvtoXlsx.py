"""
    csvファイルをxlsxファイルに変換する
    2019/01/10 Shimizu
"""
import pandas as pd
import sys

if len(sys.argv) < 2:
    print("使い方：python3 CsvtoXlsx.py xlsx化するファイル名")
    print("コマンドライン引数がありません")
    sys.exit()
args = sys.argv
csvfile_name = args[1]

#列の数を指定する
col_names = [ 'c{0:02d}'.format(i) for i in range(18) ]

# CSVファイルの読み込み
try:
    data = pd.read_csv(csvfile_name,encoding='utf-8', names=col_names )
except FileNotFoundError as e:
    print(e)
    exit()
else:
    pass

xlsmfile_name = csvfile_name.replace('csv','xlsx')
if len(sys.argv) == 3:
    data = pd.concat([data.head(1), data[1:].astype(float)])

# Excel形式で出力
data.to_excel(xlsmfile_name, encoding='utf-8',index=False,header=False)
