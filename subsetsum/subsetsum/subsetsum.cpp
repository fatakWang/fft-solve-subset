#include <iostream>
#include "fftw3.h"
#pragma comment(lib, "libfftw3-3.lib")
#include <vector>
#include<algorithm>
using namespace std;
#include <cmath>
#include<Windows.h>
fftw_complex* in_2d_fft;
fftw_complex* in_2d_ifft;
fftw_complex* out_2d_afft;
fftw_complex* out_2d_bfft;
fftw_complex* out_2d_ifft;

fftw_plan p_forward_2d_a, p_forward_2d_b, p_backward_2d;

fftw_complex* in_1d_fft, * in_1d_ifft;
fftw_complex* out_1d_afft;
fftw_complex* out_1d_bfft;
fftw_complex* out_1d_ifft;

fftw_plan p_forward_1d_a, p_forward_1d_b, p_backward_1d;

vector<int> ploy1d_fft_mult(vector<int>& acoef, vector<int>& bcoef, int u, int N)//a,b输入，c输出,u为限制大小，a_max_coef，b_max_coef表示a,b的最高次数(是真实的最高次数，所以要加1)
{
	vector<int> ccoef;
	//先根据系数构造系数表示法的数组;
	for (int i = 0; i < acoef.size(); i++)
	{
		if (acoef[i] <= u)in_1d_fft[acoef[i]][0] = 1;
	}
	fftw_execute(p_forward_1d_a);

	for (int i = 0; i < acoef.size(); i++)
	{
		if (acoef[i] <= u)in_1d_fft[acoef[i]][0] = 0;
	}
	for (int i = 0; i < bcoef.size(); i++)
	{
		if (bcoef[i] <= u)in_1d_fft[bcoef[i]][0] = 1;
	}
	fftw_execute(p_forward_1d_b);

	for (int i = 0; i < bcoef.size(); i++)
	{
		if (bcoef[i] <= u)in_1d_fft[bcoef[i]][0] = 0;
	}
	

	//已经有a，b的点值表示了接下来O（n）的时间算乘法
	for (int i = 0; i < N; i++) {  //(a+bi)(c+di) = (ac - bd) + (bc + ad)i
		in_1d_ifft[i][0] = out_1d_afft[i][0] * out_1d_bfft[i][0] - out_1d_afft[i][1] * out_1d_bfft[i][1];
		in_1d_ifft[i][1] = out_1d_afft[i][1] * out_1d_bfft[i][0] + out_1d_afft[i][0] * out_1d_bfft[i][1];
	}
	fftw_execute(p_backward_1d);
	//std::cout << "a * b poly \n";
	for (int i = 0; i < N && i <= u; i++) {
		//cout << out[i][0] / N << "\n";
		if (out_1d_ifft[i][0] / N >= 0.8) {
			ccoef.push_back(i);
		}
	}//out[i][0]表示了x^i的系数

	return ccoef;
}

void ploy2d_fft_mult(vector<vector<int> >& acoef, vector<vector<int> >& bcoef, vector<vector<int> >& ccoef, int u, int row, int col)
{
	for (int i = 0; i < acoef[0].size(); i++)
	{
		if (acoef[0][i] <= u) {
			in_2d_fft[acoef[0][i] + col * acoef[1][i]][0] = 1;
		}
	}
	fftw_execute(p_forward_2d_a);
	for (int i = 0; i < acoef[0].size(); i++)
	{
		if (acoef[0][i] <= u) {
			in_2d_fft[acoef[0][i] + col * acoef[1][i]][0] = 0;
		}
	}

	//b,再对b进行二维傅里叶变换
	for (int i = 0; i < bcoef[0].size(); i++)
	{
		if (bcoef[0][i] <= u) {
			in_2d_fft[bcoef[0][i] + col * bcoef[1][i]][0] = 1;
		}
	}

	fftw_execute(p_forward_2d_b);

	for (int i = 0; i < bcoef[0].size(); i++)
	{
		if (bcoef[0][i] <= u) {
			in_2d_fft[bcoef[0][i] + col * bcoef[1][i]][0] = 0;
		}
	}

	vector<vector<vector<double> >> bfft(row, vector<vector<double> >(col, vector<double>(2)));//三维数组，row行，col列，每个还有0，1两个范围，0表示实数，1表示虚数

	// 计算结果存储到afft中


	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < col; ++j)
		{
			in_2d_ifft[j + col * i][0] = out_2d_afft[j + col * i][0] * out_2d_bfft[j + col * i][0] - out_2d_afft[j + col * i][1] * out_2d_bfft[j + col * i][1];//(a+bi)(c+di)=ac-bd+i(ad+bc)
			in_2d_ifft[j + col * i][1] = out_2d_afft[j + col * i][0] * out_2d_bfft[j + col * i][1] + out_2d_afft[j + col * i][1] * out_2d_bfft[j + col * i][0];
		}
	}

	// ifft

	fftw_execute(p_backward_2d);
	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < col && j <= u; ++j)
		{
			if (out_2d_ifft[j + col * i][0] / (row * col) >= 0.8) {//防止浮点数的上下浮动，所以不能取1
				ccoef[0].push_back(j);
				ccoef[1].push_back(i);
			}
		}
	}




}

vector<vector<int> > all_subset_sums_withinfo(vector<int > S, int u, int row, int col)
{
	if (S.size() == 1)//表示集合仅有一个元素
	{
		vector<vector<int> > re(2, vector<int>(2));
		re[0][0] = re[1][0] = 0;
		re[0][1] = S[0];
		re[1][1] = 1;
		return re;
	}
	else if (S.size() == 0) {//表示为空集
		vector<vector<int> > re(2, vector<int>(1));
		re[0][0] = re[1][0] = 0;
		return re;
	}
	else {
		int size = S.size();
		vector<int> T(S.begin(), S.begin() + size / 2);
		vector<int> S_div_T(S.begin() + size / 2, S.end());
		vector<vector<int> > T_tem = all_subset_sums_withinfo(T, u, row, col);
		vector<vector<int> > S_div_T_tem = all_subset_sums_withinfo(S_div_T, u, row, col);
		vector<vector<int> > re(2);
		ploy2d_fft_mult(T_tem, S_div_T_tem, re, u, row, col);
		return re;
	}
}

vector<int > all_subset_sums(vector<int > S, int u)
{
	int n = S.size();
	int b = sqrt(log(n) * n);
	if (b == 0) {
		vector<int > re = { 0 };
		re = ploy1d_fft_mult(re, S, u, 0);
		return re;
	}
	int l = 0;
	vector<vector<int> > Q_l(b);
	vector<vector<int> > R_l(b);
	for (int i = 0; i < n; i++)
	{
		l = S[i] % b;//l表示除b的余数
		Q_l[l].push_back(S[i] / b);//在第l行加入与b的除数
	}//构造Q_l
	int row = n + 1, col = (u / b) * 2 + 1;

	in_2d_fft = (fftw_complex*)fftw_malloc(row * col * sizeof(fftw_complex));
	in_2d_ifft = (fftw_complex*)fftw_malloc(row * col * sizeof(fftw_complex));
	out_2d_afft = (fftw_complex*)fftw_malloc(row * col * sizeof(fftw_complex));
	out_2d_bfft = (fftw_complex*)fftw_malloc(row * col * sizeof(fftw_complex));
	out_2d_ifft = (fftw_complex*)fftw_malloc(row * col * sizeof(fftw_complex));

	p_forward_2d_a = fftw_plan_dft_2d(row, col, in_2d_fft, out_2d_afft, FFTW_FORWARD, FFTW_MEASURE);
	p_forward_2d_b = fftw_plan_dft_2d(row, col, in_2d_fft, out_2d_bfft, FFTW_FORWARD, FFTW_MEASURE);
	p_backward_2d = fftw_plan_dft_2d(row, col, in_2d_ifft, out_2d_ifft, FFTW_BACKWARD, FFTW_MEASURE);

	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < col; ++j)
		{
			in_2d_fft[j + col * i][0] = 0;
			in_2d_ifft[j + col * i][0] = 0;
			in_2d_fft[j + col * i][1] = 0;
			in_2d_ifft[j + col * i][1] = 0;
		}
	}



	for (l = 0; l < b; l++)
	{
		vector<vector<int> > S_Q_l = all_subset_sums_withinfo(Q_l[l], u / b, row, col);
		int x;
		for (int i = 0; i < S_Q_l[0].size(); i++)
		{
			x = S_Q_l[0][i] * b + S_Q_l[1][i] * l;
			R_l[l].push_back(x);
		}
	}//构造R_l
	vector<int > re = R_l[0];
	int N = 2 * u + 1;
	in_1d_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	in_1d_ifft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out_1d_afft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out_1d_bfft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	out_1d_ifft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
	for (int i = 0; i < N; i++) {
		in_1d_fft[i][0] = 0;
		in_1d_fft[i][1] = 0;
		in_1d_ifft[i][0] = 0;
		in_1d_ifft[i][1] = 0;
	}
	p_forward_1d_a = fftw_plan_dft_1d(N, in_1d_fft, out_1d_afft, FFTW_FORWARD, FFTW_MEASURE);
	p_forward_1d_b = fftw_plan_dft_1d(N, in_1d_fft, out_1d_bfft, FFTW_FORWARD, FFTW_MEASURE);
	p_backward_1d = fftw_plan_dft_1d(N, in_1d_ifft, out_1d_ifft, FFTW_BACKWARD, FFTW_MEASURE);
	for (l = 0; l < b - 1; l++)
	{
		re = ploy1d_fft_mult(re, R_l[l + 1], u, N);
	}
	fftw_destroy_plan(p_forward_2d_a);
	fftw_destroy_plan(p_forward_2d_b);
	fftw_destroy_plan(p_backward_2d);
	fftw_destroy_plan(p_forward_1d_a);
	fftw_destroy_plan(p_forward_1d_b);
	fftw_destroy_plan(p_backward_1d);
	return re;
}

int main() {

	LARGE_INTEGER t1, t2, tc;
	QueryPerformanceFrequency(&tc);
	QueryPerformanceCounter(&t1);
	vector<int> S = { 518533,1037066,2074132,1648264,796528,1593056, 686112,1372224,244448,488896,977792,1955584,1411168,322336,644672,1289344,78688,157376,314752,629504,1259008 };
	vector<int > re = all_subset_sums(S, 20000000);//行取值为0，1，0表示x的次方，1表示y的次方。列取值为点的个数,在这里x是子集的和，y是对应的子集的元素个数
	QueryPerformanceCounter(&t2);
	double time = (double)(t2.QuadPart - t1.QuadPart) / (double)tc.QuadPart;
	cout << "time = " << time << endl;  //输出时间（单位：ｓ）
	for (int i = 0; i < re.size(); i++)
	{
		cout << re[i] << endl;
	}
	//text();
	int i;
	cin >> i;

}