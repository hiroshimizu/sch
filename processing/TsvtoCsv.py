"""
    tab & space を , に置換してcsvファイルを新たに出力する
    2019/01/09 Shimizu
"""
import re, sys

if len(sys.argv) < 2:
    print("使い方：python3 TsvtoCsv.py csv化するファイル名")
    print("コマンドライン引数がありません")
    sys.exit()

args = sys.argv
beginning_tab_and_space = re.compile(r'^[\t, ]+')
tab_and_space = re.compile(r'[\t, ]+')
tsvfile_name = args[1]

try:
    tsvfile = open(tsvfile_name,"r", encoding="utf-8")
except FileNotFoundError as e:
    print(e)
    exit()
else:
    pass

csvfile_name = tsvfile_name.replace('.','_') + '.csv'
csvfile = open(csvfile_name, "w", encoding="utf-8")
r = tsvfile.readlines()
for row in r:
    csvfile.write(tab_and_space.sub(',',beginning_tab_and_space.sub('',row)))
    
tsvfile.close()
csvfile.close()
