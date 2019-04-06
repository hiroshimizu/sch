#!/bin/sh
#*************************************************************************
# execute.sh �f�[�^�t�@�C�����G�N�Z���Ɏ����Ŏ�荞�ނ��߂̃v���O�����Q��
#				 execute�i���s����j�V�F���X�N���v�g shimizu 2019/01/12
#*************************************************************************

#*************************************************************************
#���e
# 1.�X�N���v�g�Q�����s����f�B���N�g���i���j�ֈ�U�R�s�[����(4.�ŏ���)
# 2.�f�[�^�t�@�C��(tsv)��csv���o�R����xlsx�t�@�C���֕ϊ�
# 3.xlsx�t�@�C�����e���v���[�g�֏�������
# 4.����Ȃ��t�@�C��������
#*************************************************************************

########################################################################1.
cp copy_sheet.py ../
cp copy_text_to_xlsx.py ../
cp CsvtoXlsx.py ../
cp TsvtoCsv.py ../
cd ../
echo 'Coping neccessary files has been completed.'
#�����̃R�����g�L���ŗv�f�����Ȃ��悤�ɃX�y�[�X������
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
