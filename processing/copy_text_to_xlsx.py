"""
    テキストファイルの内容をxlsxファイルのシートにコビーする
    使い方：python3 copy_text_to_xlsx.py txtfile wb sheet
    2019/01/10 Shimizu
"""
import openpyxl, sys

if len(sys.argv) != 4:
    print("使い方：python3 copy_text_to_xlsx.py txtfile wb sheet")
    print("コマンドライン引数は3つ:txtfile,book,sheetの順")
    sys.exit()

args = sys.argv
txt_name = args[1]
wb_name = args[2]
sheet_name = args[3]

txtfile = open(txt_name)

try:
    wb = openpyxl.load_workbook(wb_name)
except FileNotFoundError as e:
    print(e)
    exit()
else:
    pass

try:
    sheet = wb.get_sheet_by_name(sheet_name)
except KeyError as e:
    print(e)
    exit()
else:
    pass
i=1
for line in txtfile:
    sheet.cell(row=i,column=1).value = line
    i += 1

wb.save(wb_name)
txtfile.close()
