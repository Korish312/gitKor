#include <stdio.h>
#define N 12

double detA(double matA[N][N], int n);
void cofactorA(double matA[N][N], int n, double ans[N][N]);
void multiply(double matA[N][N], double a, double astar[N][N],int n);
int rankA(double matA[N][N], int n);

int main()
{
	double matA[N][N];//matA矩阵
	double astar[N][N];//A*
	/*,需要注意的一点是，这里最好用double而不要用int，因为一个问题是矩阵的元素可以是小数
另一个是我们在接下来求逆矩阵和两矩阵相乘验证是否为单位阵的时候，用int可能会造成truncate而导致结果不对*/
	int i, j;
	int n;//矩阵阶数n

	printf("please enter the order of the matrix:");
	scanf_s("%d", &n);
	while (n < 1|| n > N)//判断n是否有效
	{
		printf("Error!!");
		return 0;
	}
	while (n>=1 && n <= N)
	{
		printf("please enter the elements of the matrix:");
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				scanf_s("%lf", &matA[i][j]);//输入矩阵
			}
		}

		double a = detA(matA, n);
		if (a == 0)//利用detA检验A是否可逆，若detA=0，矩阵不可逆
		{
			printf("There's no inverse matrix!\n");
			printf("是否求该矩阵的秩？\n");
			printf("若是，请输入1，若否，输入其他值:");
			int c = 0;
			scanf_s("%d", &c);
			if (c == 1)
			{
				rankA(matA, n);
				return 0;
			}
			else
			{
				return 0;
			}
			return 0;
		}
		else
		{
			printf("矩阵的秩为%d\n", n);
			printf("逆矩阵为\n");
			cofactorA(matA, n, astar);
			for (i = 0; i < n; i++)
			{
				for (j = 0; j < n; j++)
				{
					printf("%.4lf ",astar[i][j] / a);//利用公式A^-1=(CofactorA/det A)求逆矩阵
				}
				printf("\n");
			}

		printf("是否需要验证逆矩阵与原矩阵的积为单位矩阵\n");//这一块用于验证两矩阵相乘是否为单位矩阵
		printf("若是，请输入1,若否，输入其他值：");
		int b = 0;
		scanf_s("%d", &b);
		if (b == 1)
		{
			multiply(matA, a, astar,n);
		}
		return 0;
		}
		printf("\n");
	}
	return 0;
}

double detA(double matA[N][N], int n)//按第一行展开计算A的行列式
{
	if (n == 1)//n=1时，行列式的值等于数字的值
	{
		return matA[0][0];
	}
	double ans = 0;//用ans来记录各项的加减
	double temp[N][N];//定义临时的矩阵来记录cofactor
	int i, j, k;
	for (k = 0; k < n; k++)//对第一行的第k个元素，计算其cofactor
	{
		for (i = 0; i < n - 1; i++)
		{
			for (j = 0; j < n - 1; j++)
			{
				temp[i][j] = matA[i + 1][(j >= k) ? j + 1 : j];
				/*关于这里的注释是，由于我们是用第一行展开计算的行列式，所以cofactor的第i行就相当于A的第i+1行
				而cofactor的列则是扒掉第k列之后（因为计算的是第k个元素的cofactor），若j小于k,则还是原来A中的列数不变，若大于等于k，则需要加一*/
			}
		}
		double t = detA(temp, n - 1);//这里我们用int t来表示我们刚刚算出来的cofactor矩阵的行列式（其实这里函数名用detA似乎并不是特别合适。。。）
		if (k % 2 == 0)//注意到求行列式的公式里涉及正负号的问题。我们在这里用matA的第一行的第k个元素乘以其det cof
		{
			ans += matA[0][k] * t;
		}
		else
		{
			ans -= matA[0][k] * t;
		}
	}
	return ans;
}

void cofactorA(double matA[N][N], int n, double ans[N][N])//计算每一行每一列的每个元素所对应的余子式，组成A*
{
	if (n == 1)//对n=1阶矩阵(其实就是一个数），其cofactor等于1
	{
		ans[0][0] = 1;
		return;
	}
	int i, j, k, t;//上面我们是计算第一行的cofactor并计算其det使用，此次我们要计算每一行每一列的，所以需要四个变量
	double temp[N][N];//定义一个临时的矩阵，用于存储计算出的cofactor的结果
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			for (k = 0; k < n - 1; k++)
			{
				for (t = 0; t < n - 1; t++)
				{
					temp[k][t] = matA[k >= i ? k + 1 : k][t >= j ? t + 1 : t];
					/*跟上面detA函数里面基本上是一样的操作，不同的地方是这里加入了列的问题*/
				}
			}


			ans[j][i] = detA(temp, n - 1);//在这里我们保存了计算出的cofactor matrix，需要注意的一点是这里是ans[j][i]而不是ans[i][j]，因为这里要转置
			if ((i + j) % 2 == 1)//同上面一样的正负号判据，但此处需同时考虑i和j二者的值
			{
				ans[j][i] = -ans[j][i];
			}
		}
	}
}

void multiply(double matA[N][N], double a, double astar[N][N],int n)
{
	int i, j;
	double C[N][N];
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			long double sum = 0;
			for (int m = 0; m < n; m++)
			{
				sum = sum + matA[i][m] * astar[m][j] / a;
			}
			C[i][j] = sum;//矩阵乘法算出C矩阵的ij元素
		}
	}
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf("%.1lf ", C[i][j]);
		}
		printf("\n");
	}
}

int rankA(double matA[N][N], int n)
{
	double temp[N][N];
	int i, j, r;
	int a = n;
	int b = 0;

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			temp[i][j] = matA[i][j];
		}
	}
	while (detA(temp, n) == 0)
	{
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				temp[i][j] = temp[i+1][j+1];
			}
		}
		n--;
		b++;
	}
	r = a - b;
	printf("矩阵的秩为%d", r);
	return r;
} 
