#pragma once  // ヘッダファイルにはこれを書く

#ifndef __SCHRO_CONSTANT__
#define __SCHRO_CONSTANT__
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<limits>
#include<map>
#include<cstdio>

const REAL pi = 3.1415926535897932384626433832795028841971693993751;


// 文字列の空白を取り除く
std::string trim(std::string str)
{	std::stringstream ss(str);
	std::string s;
	ss >> s;
	return s;
}

#define READ_VALUE(NAME,LIST) \
	NAME=read_value(NAME,#NAME,LIST)
template<class T>
T read_value(
	T& var, std::string name,
	std::map<std::string,std::string>& param_list)
{	typedef std::map<std::string,std::string> list_type;

	// 読み込んだパラメータリストにあるか探す
	const list_type::iterator itr=param_list.find(name);
	if(itr!=param_list.end()) {
		// 適切な型に変換して代入
		std::string str = itr->second;
		std::stringstream ss(str);
		ss >> var;
	}
	else {
		// リストになければNaNにしておく(整数の場合は0)
		var = std::numeric_limits<T>::quiet_NaN();
		std::cout << "Error: Parameter is not found in file" << std::endl;
	}
	return var;
}

std::map<std::string,std::string>
read_parameter_list(std::string filename)
{	using std::string;
	std::ifstream ifs(filename.c_str());
	std::map<string,string> parameter_list;
	while(ifs) {
		// 一行読み込み
		string line;
		std::getline(ifs, line);

		// # 以下はコメント扱いとして消去
		string::size_type pos=line.find("#");
		if(pos!=string::npos) {
			line.erase(line.begin()+pos,line.end());
		}

		// "="で分割
		std::stringstream ss(line);
		string first, second;
		std::getline(ss, first,'=');
		std::getline(ss,second,'=');

		// 単語前後のスペースを除去
		first = trim(first);
		second= trim(second);

		// 取得結果を保存
		if(!first.empty() && !second.empty()) {
			parameter_list[first] = second;
		}
	}
	return parameter_list;
}

bool read_parameter(std::string filename,
		int& N_rec, int& N_rec_rho, int& N_step, REAL& hx, REAL& hy, REAL& Lx, REAL& Ly,
		REAL& omg, REAL& v0x, REAL& v0y, REAL& Bz, REAL& m, REAL& q, REAL& x0, REAL& y0,
		REAL& vSIZE, REAL& v0_, REAL& Bz_, REAL& mp_, REAL& e_, REAL& h_bar_)
{	// 定数読み込み
	std::map<std::string,std::string>
		parameter_list = read_parameter_list(filename);
	if(parameter_list.empty())return false;

	// (READ_VALUEは宣言していない変数を使うとエラーとなる)
	READ_VALUE(N_rec	 ,parameter_list);
	READ_VALUE(N_rec_rho ,parameter_list);
	READ_VALUE(N_step	 ,parameter_list);
	READ_VALUE(hx      	 ,parameter_list);
	READ_VALUE(hy      	 ,parameter_list);
	READ_VALUE(Lx      	 ,parameter_list);
	READ_VALUE(Ly      	 ,parameter_list);
	READ_VALUE(omg    	 ,parameter_list);
	READ_VALUE(v0x    	 ,parameter_list);
	READ_VALUE(v0y    	 ,parameter_list);
	READ_VALUE(Bz		 ,parameter_list);
	READ_VALUE(m		 ,parameter_list);
	READ_VALUE(q		 ,parameter_list);
	READ_VALUE(x0   	 ,parameter_list);
	READ_VALUE(y0   	 ,parameter_list);
	READ_VALUE(vSIZE  	 ,parameter_list);
	READ_VALUE(v0_  	 ,parameter_list);
	READ_VALUE(Bz_    	 ,parameter_list);
	READ_VALUE(mp_  	 ,parameter_list);
	READ_VALUE(e_     	 ,parameter_list);
	READ_VALUE(h_bar_ 	 ,parameter_list);

	return true;
}
#endif
