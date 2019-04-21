#**********************************************************************************************
# set conditions.py
# conditions.pyを書き換えるプログラム
# 原点を中心に円運動させるためサイクロトロン半径 x0=mv/qB に注意
#{'Case':'001', 'v0x':'0', 'v0y':'1', 'B':'1', 'm':'1', 'q':'1', 'x0':'-1','y0':'0', 'h':'1', 'E':'1e-4'}
#**********************************************************************************************
import pprint

file_obj = open('conditions.py', 'w')
Conditions = []

y0 = str(0)
v0x = str(0)
case_number = 0
#for count_v0x in range(1, 3):
#	v0x = str(count_v0x)
for count_v0y in range(1, 3):
	v0y = str(count_v0y)
	for count_B in range(1, 3):
		B = str(count_B)
		for count_m in range(1, 3):
			m = str(count_m)
			for count_q in range(1, 3):
				q = str(count_q)
				x0_num = count_v0y * count_m / count_B / count_q
				x0 = str(x0_num)
#				for count_x0 in range(3):
#					x0 = str(count_x0)
#					for count_y0 in range(3):
#						y0 = str(count_y0)
				for count_h in range(1, 7):
					h = str(count_h)
					for count_E in range(1,3):
						E = '1e-' + str(count_E)
						case_number += 1
						condition = {'Case':str(case_number), 'v0x':v0x, 'v0y':v0y, 'B':B, 'm':m, 'q':q, 'x0':x0,'y0':y0, 'h':h, 'E':E}
						Conditions.append(condition)

file_obj.write('Conditions = ' + pprint.pformat(Conditions))
file_obj.close()
