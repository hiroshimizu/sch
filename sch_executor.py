#***********************************************************************************************
#									sch_executor.py
# 条件毎のディレクトリを作成し、sch.inp potential.cuをそれぞれ書き換え、すべて実行する
# Dataディレクトリにすべてのデータ、resultディレクトリに必要最小限のデータが入る
#使い方：conditions.pyに条件を書いた後にpython3 sch_executor.py を実行する
#***********************************************************************************************
import re, sys, os, conditions, logging, shutil, subprocess
logging.basicConfig(level=logging.DEBUG,format=' %(asctime)s - %(levelname)s - %(message)s')
logging.disable(logging.CRITICAL)
#**********************************************************************************************
#											関数get_filename(cdtn)
#         引数：辞書型の条件	戻り値：Case別フォルダ(Caseフォルダ)名の文字列
#**********************************************************************************************
def get_filename(cdtn):
	filename = 'Case' + cdtn['Case'] + 'u' + cdtn['v0x'] + 'v' + cdtn['v0y'] + 'B' + cdtn['B'] + \
			'm' + cdtn['m'] + 'q' + cdtn['q'] + 'x' + cdtn['x0'] + 'y' + cdtn['y0'] + \
			'h' + cdtn['h'] + 'E' + cdtn['E']
	return filename
#****************************************************************************
#											ufm_judge
# 						ufmを計算するかどうか判定する関数
#				E以外全て同じ条件のCaseがある場合、ufmはその結果を使えばよいので計算が省ける
#              この関数でその結果もコピーもする
# ufmを実行する場合はufm_judge()はTureを、実行しない場合はFalseを戻り値にとる
#****************************************************************************

def ufm_judge(i):
	boolian = True
	for k in range(i):
		if conditions.Conditions[k]["v0x"]==conditions.Conditions[i]["v0x"] \
			and conditions.Conditions[k]["v0y"]==conditions.Conditions[i]["v0y"] \
			and conditions.Conditions[k]["B"]==conditions.Conditions[i]["B"] \
			and conditions.Conditions[k]["m"]==conditions.Conditions[i]["m"] \
			and conditions.Conditions[k]["q"]==conditions.Conditions[i]["q"] \
			and conditions.Conditions[k]["x0"]==conditions.Conditions[i]["x0"] \
			and conditions.Conditions[k]["y0"]==conditions.Conditions[i]["y0"] \
			and conditions.Conditions[k]["h"]==conditions.Conditions[i]["h"]:
				print("{}と{}はufmが同じです".format(conditions.Conditions[i]["Case"], conditions.Conditions[k]["Case"]))
				os.makedirs('./Data/' + get_filename(conditions.Conditions[i]) + "_ufm")
				os.system('touch '+ get_filename(conditions.Conditions[i]) + '/equal' + get_filename(conditions.Conditions[k]))
				shutil.copy('./Data/' + get_filename(conditions.Conditions[k]) + '_ufm/000.ufm', './Data/' + get_filename(conditions.Conditions[i]) + '_ufm'))
				boolian = False
	return boolian

#**********************************************************************************************
#											関数get_fileObj(filepath)
#         引数：ファイルのパス（文字列）	戻り値:ファイルオブジェクト
#			 機能:try-except構文を備えたopen関数
#**********************************************************************************************
def get_fileObj(filepath):
	try:
		fileObj = open(filepath,'r',encoding='utf-8')
	except Exception as err:
		print(filepath + 'が開けませんでした:' + str(err))
		exit()
	return fileObj
#**********************************************************************************************
#											関数change_condition
#         引数：辞書型の条件	機能:Caseフォルダ内のpotential.cuとsch.inpを書き換える
#**********************************************************************************************
def change_condition(condition_dict):
	file_name = 'Data/' + get_filename(condition_dict)
	#方針：メモリ上にpotential.cuとsch.inpのテキストをすべて乗せ、regex.sub()で書き換える
	#		 一時ファイルにその書き換えたデータを出力する
	potential_cu = get_fileObj(file_name + '/potential.cu')
	sch_inp = get_fileObj(file_name + '/sch.inp')
	temp_potential_cu = open(file_name + '/temp_potential.cu','w',encoding='utf-8')
	temp_sch_inp = open(file_name + '/temp_sch.inp','w',encoding='utf-8')
	#potential.cuとsch.inpをメモリ上に乗せる
	potential_lines = potential_cu.readlines()
	sch_lines = sch_inp.readlines()
	target_regex = re.compile(r'INP_NUM_')
	#メモリ上のpotential.cuを書き換える
	for i in range(0,len(potential_lines)):
		mo = target_regex.search(potential_lines[i])
		if mo != None:
			potential_lines[i] =  target_regex.sub(condition_dict['E'],potential_lines[i])
			break
	#メモリ上のsch.inpを書き換える
	for parameter in condition_dict.keys():
		if parameter.startswith('Case'):
			continue
		parameter_regex = re.compile('{}'.format(parameter))
		for i in range(0,len(sch_lines)):
			mo = parameter_regex.search(sch_lines[i][0:5])
			if mo != None:
				if parameter == 'h':
					sch_lines[i] =  target_regex.sub((str(float(condition_dict[parameter])*1.054) + 'e-34').rjust(11),sch_lines[i])
				else:
					sch_lines[i] =  target_regex.sub(condition_dict[parameter].rjust(8),sch_lines[i])
	#メモリ上のデータを一時ファイルに出力する
	for i in range(0,len(potential_lines)):
		temp_potential_cu.write(potential_lines[i])
		#logging.debug('i={},{}'.format(i,potential_lines[i]))
	for i in range(0,len(sch_lines)):
		temp_sch_inp.write(sch_lines[i])
		#logging.debug('i={},{}'.format(i,sch_lines[i]))
	potential_cu.close()
	temp_potential_cu.close()
	sch_inp.close()
	temp_sch_inp.close()
	shutil.move(file_name + '/temp_potential.cu',
			file_name + '/potential.cu')
	shutil.move(file_name + '/temp_sch.inp',
			file_name + '/sch.inp')
#****************************************************************************
#										make
#****************************************************************************
def command_make(file_path):
	os.chdir(file_path)
	os.makedirs("./OUT")
	os.makedirs("./out")
	os.makedirs("./SAVE")
	cuda_capability = subprocess.check_output("deviceQuery | grep CUDA\ Capability | cut -d ' ' -f 11", shell=True).decode()[:3]
	if cuda_capability == "2.0":
		print("cuda_capability == 2.0")
		shutil.copy("Makefile.cc20", "Makefile")
	elif cuda_capability == "5.2" or cuda_capability == "6.1":
		print("cuda_capability == 5.2 OR 6.1")
		shutil.copy("Makefile.cc61", "Makefile")
	elif cuda_capability == "7.5":
		print("cuda_capability == 7.5")
		shutil.copy("Makefile.cc75", "Makefile")
	else:
		print("エラー:deviceQueryでCUDA Capabilityを確認し、適切なMakefileを用意してください。")
		exit(1)
	os.system('make')
	if file_path.endswith('ufm'):
		os.system('./sch_ufm > aaa.ufm')
	else:
		os.system('./sch > aaa.out')
		os.chdir('./PLT')
		os.system('gnuplot toEPS_multi3.plt')
		os.chdir('../')
	os.chdir('../../')

#****************************************************************************
#										Main
#****************************************************************************
if "Data" not in os.listdir("."):
	os.makedirs('./Data')
if "result" not in os.listdir("."):
	os.makedirs('./result')

for i in range(len(conditions.Conditions)):
	file_name = 'Data/' + get_filename(conditions.Conditions[i])
	if os.path.exists(file_name):
		continue
	print('making ' + file_name +'...')
	shutil.copytree('Source',file_name)
	change_condition(conditions.Conditions[i])
	if(ufm_judge(i)):
		shutil.copytree(file_name, file_name + '_ufm')
		command_make(file_name + '_ufm')
	command_make(file_name)
	file_name2 = 'result/' + get_filename(conditions.Conditions[i])
	os.makedirs(file_name2)
	shutil.copy(file_name + '/000.dat', file_name2)
	shutil.copy(file_name + '_ufm/000.ufm', file_name2)
	shutil.copy(file_name + '/PLT/locus_info.eps', file_name2)
	shutil.copy(file_name + '/potential.cu', file_name2)
	shutil.copy(file_name + '/sch.inp', file_name2)
