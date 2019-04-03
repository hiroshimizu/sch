#include <stdio.h>
#include <string.h>

/***** 波動関数をbinary形式で書き出す関数 *************************/
/*	binary形式でPsiを書き出す                                      */
/*	Input : 	const int k, 書き出すファイルの名前にするループ回数k。*/
/* Output:	CMPLX* Psi, 波動関数の配列                            */
/* 			REAL E0, その波動関数を得たときの初期エネルギー       */
/*		呼び出し元で                                                */
/*		if(k==300) save_wavefunction(k, Psi, E0);                   */
/*		などと使うことを想定している。                              */
/******************************************************************/
void save_wavefunction(const int NG, const int k, CMPLX* Psi, REAL E0)
{	printf("Now Saving ...");
	char str_Psi[20]; sprintf(str_Psi,"SAVE/%08d.bin",k);
	char str_EN0[] = "SAVE/EN0.bin";
	//printf("check str_EN0 = %s\n", str_EN0);
	FILE *fPsi, *fEN0;
	if (( fPsi = fopen(str_Psi, "wb") ) == NULL )
	{	printf(" file open error!!__%s", str_Psi);	}
	else 
	{	fwrite(Psi, sizeof(CMPLX), NG, fPsi); fclose(fPsi);}
	if (( fEN0 = fopen(str_EN0, "wb") ) == NULL )
	{ 	printf(" file open error!!__%s", str_EN0);	}
	else
	{	fwrite(&E0, sizeof(REAL),   1, fEN0); fclose(fEN0); }
	//	バイナリファイルを書き出したらkのループを抜けて終了処理を行う
	printf(" done!\n");
}

/**** 波動関数を読み込む関数 *********************************/
/*	ディレクトリSAVEに格納されているに"*.bin"の中でもっとも   */
/*	値が大きいものを波動関数の配列に読み込む                  */
/*		Input :	int N_t_step, 計算全体のforループの周回数     */
/*					int NG, 読み込む波動関数の全節点数            */
/*					CMPLX*, Psi 読み込む波動関数                  */
/*					REAL* E0, 読み込む波動関数の初期エネルギー    */
/*		Output:	int (戻り値), loadできたかどうか、0:偽、1:真  */
/*					int* k_start, restartするときのkの値          */
/*	波動関数が発見されればPsiとE0に値が読み込まれる           */
/*	発見されなければ、Psi, E0,らに対して何もしない            */
/*************************************************************/
int load_wavefunction(const int N_t_step, const int NG, int* k_start, CMPLX* Psi, REAL* E0)
{ /* SAVE内の"000~~~~~.bin"を検索して、見つかったらそこからrestartする */
	FILE *fPsi, *fEN0;
	char str_EN0[] = "SAVE/EN0.bin";
	char str_Psi[20];
	printf("Searcing for SAVE/EN0.bin\n");
	if ( (fEN0 = fopen(str_EN0, "rb")) == NULL )
	{	printf("Initial Energy is NOT found.\n"); 
		return 0; // False
	}
	for (int k = N_t_step; k >= 0; k--)	// kの値が大きい方から照合していく
	{	sprintf(str_Psi, "SAVE/%08d.bin", k);
		printf("Searching for %s\r", str_Psi);
		if ( (fPsi = fopen(str_Psi, "rb")) != NULL )
		{	//	fopenの戻り値がNULLでなかったならば、
			//	すなわちfileが見つかったならば
			printf("Savedata is found: %s\n", str_Psi);
			fread(Psi, sizeof(CMPLX), NG, fPsi);
			fread(E0,  sizeof(REAL),   1, fEN0);
			fclose(fPsi);fclose(fEN0);
			printf(" done!\n");
			*k_start = k + 1;	// savedataの次の周回から始める
			return 1; // Ture
		}
	}	
	printf("00~~~~~~.bin file is NOT found.\n");
	return 0; // False
}
