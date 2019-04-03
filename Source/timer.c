/******	msの単位で経過時間を測る関数(開始用)	****/
/*	入力:	StopWatchInterface型のポインタポインタ	   */
/*	出力:	なし                                      */
/*			単なる準備用。使うときはメインの方で      */
/*			 StopWatchInterface* timer;               */
/*			などとタイマー用の変数を作って            */
/*			 SetTimer(&timer);                        */
/*			と & をつけて引数に書く                   */
/***************************************************/
void SetTimer(StopWatchInterface** TIMER)
{	
	sdkCreateTimer(TIMER);
	sdkStartTimer(TIMER); 
}

/******	msの単位で経過時間を測る関数(終了用)********/
/*	入力:	StopWatchInterface型のポインタポインタ     */
/*	出力:	(戻り値)double型でSetTimerが使われてからの */
/*			経過時間を得る                             */
/*			 時間取得兼タイマー終了用。使うときは      */
/*			 double time = EndTimer(&timer);           */
/*			と & をつけて引数に書く                    */
/****************************************************/
double EndTimer(StopWatchInterface** TIMER)
{	double tmp = -1;
	sdkStopTimer(TIMER);
	tmp = sdkGetTimerValue(TIMER);
	sdkDeleteTimer(TIMER);
	return tmp;
}
