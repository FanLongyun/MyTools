#pragma once

#include <math.h>
#include <exception>
#include <stdexcept>
#include <vector>
#include <algorithm>

#define PI	3.1415926535897932384626
#include "coordcnv.h"

#ifndef DLL_API
#define DECLDIR __declspec(dllimport)
#else
#define DECLDIR __declspec(dllexport)
#endif

#define DEG2RAD(deg) deg * PI / 180.0
#define RAD2DEG(rad) rad * 180.0 / PI

struct POSITIONSTRUCT
{
	double latitude;
	double longitude;
	double altitude;
	float fHeading;
	float fPitch;
	float fRoll;
};

// 计算目标与本机相对航向
DECLDIR double Angle(const double &destLat, const double &destLon, const POSITIONSTRUCT *P3DData);

// 线性变换算法
DECLDIR double LinearTransform(const double& input, const double& InStart, const double& InEnd, const double& OutStart, const double& OutEnd);

// 矩阵类
template <typename T>
class Matrix
{
public:
	Matrix();
	Matrix(int row, int column, T num);		// row:行数 column:列数 num:元素初始值
	Matrix(const Matrix<T>&);
	~Matrix();
	Matrix<T> operator *(const Matrix<T>&) const;
	Matrix<T> operator +(const Matrix<T>&) const;
	Matrix<T> operator -(const Matrix<T>&) const;
	Matrix<T> operator =(const Matrix<T>&);
	template <typename G>
	Matrix<G> operator / (const G& scalar) const;
	Matrix<T> Transposition();					// 转置矩阵
	T Determinant() const; 						// 行列式值
	Matrix<T> Adjoint() const;					// 标准伴随矩阵
	template <typename G>
	Matrix<G> Inverse(G&) const;					// 逆矩阵
	T** element;
	int row;
	int column;
};
// 生成矩阵
template <typename T>
Matrix<T> MakeMatrix(int row, int column, T num);

// 计算方向余弦矩阵
template <typename T>
Matrix<T> SetDCM(float psi, float theta, float phi);	// 单位:角度  psi: -180 --- 180 右转为正  theta: -90 --- 90 向上为正

// 矩阵转置
template <typename T>
Matrix<T> MatrixTranspose(Matrix<T>& mat);

template <typename T>
Matrix<T>::Matrix()
{
	row = 3;
	column = 3;
	element = new T*[row];
	for (int i = 0; i < row; i++)
		element[i] = new T[column];
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < column; j++)
			element[i][j] = 0;
	}

}

template <typename T>
Matrix<T>::Matrix(int row, int column, T num)
{
	this->row = row;
	this->column = column;
	this->element = new T*[this->row];
	for (int i = 0; i < row; i++)
		this->element[i] = new T[this->column];
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < this->column; j++)
			this->element[i][j] = num;
	}
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T>& mat)
{
	this->row = mat.row;
	this->column = mat.column;
	this->element = new T*[mat.row];
	for (int i = 0; i < mat.row; i++)
		this->element[i] = new T[mat.column];
	for (int i = 0; i < mat.row; i++)
	{
		for (int j = 0; j < column; j++)
			this->element[i][j] = mat.element[i][j];
	}
}

template <typename T>
Matrix<T>::~Matrix()
{
	for (int i = 0; i < row; i++)
		delete[]element[i];
	delete[]element;
	row = 0;
	column = 0;
}


template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& mat) const
{
	try{
		if (this->column != mat.row)
			throw std::runtime_error("row or column doesn't suit.");
	}
	catch (std::runtime_error err){
		return *this;
	}

	Matrix<T> result(this->row, mat.column, 0);
	for (int i = 0; i < this->row; i++)
	{
		for (int j = 0; j < mat.column; j++)
		{
			for (int k = 0; k < this->column; k++)
				result.element[i][j] += this->element[i][k] * mat.element[k][j];
		}
	}
	return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& mat) const
{
	try{
		if (this->row != mat.row || this->column != mat.column)
			throw std::runtime_error("row or column doesn't suit.");
	}
	catch (std::runtime_error err){
		return *this;
	}

	Matrix<T> result(mat.row, mat.column, 0);
	for (int i = 0; i < mat.row; i++)
	{
		for (int j = 0; j < mat.column; j++)
			result[i][j] = this->element[i][j] + mat.element[i][j];
	}
	return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& mat) const
{
	try{
		if (this->row != mat.row || this->column != mat.column)
			throw std::runtime_error("row or column doesn't suit.");
	}
	catch (std::runtime_error err){
		return *this;
	}

	Matrix<T> result(mat.row, mat.column, 0);
	for (int i = 0; i < mat.row; i++)
	{
		for (int j = 0; j < mat.column; j++)
			result[i][j] = this->element[i][j] - mat.element[i][j];
	}
	return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator=(const Matrix<T>&mat)
{
	try{
		if (this->row != mat.row || this->column != mat.column)
			throw std::runtime_error("row or column doesn't suit.");
	}
	catch (std::runtime_error err){
		return *this;
	}

	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < column; j++)
			this->element[i][j] = mat.element[i][j];
	}
	return *this;
}

template <typename T>
template <typename G>
Matrix<G> Matrix<T>::operator / (const G& scalar) const
{
	Matrix<G> errMat(this->row, this->column, -1.0);
	try{
		if (scalar == 0)
			throw std::runtime_error("scalar can't be 0.");
	}
	catch (std::runtime_error err)
	{
		return errMat;
	}

	Matrix<G> result(this->row, this->column, 0);
	for (int i = 0; i < this->row; i++)
	{
		for (int j = 0; j < this->column; j++)
			result.element[i][j] = static_cast<G>(this->element[i][j] / scalar);
	}

	return result;
}

template <typename T>
Matrix<T> Matrix<T>::Transposition()
{
	Matrix<T> matTemp = *this;
	for (int i = 0; i < row; i++)
		delete[]element[i];
	delete[]element;

	row = matTemp.column;
	column = matTemp.row;
	element = new T*[row];
	for (int i = 0; i < row; i++)
		element[i] = new T[column];
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < column; j++)
			element[i][j] = matTemp.element[j][i];
	}

	return *this;
}

template <typename T>
T Matrix<T>::Determinant() const
{
	try{
		if (this->row != 3 || this->column != 3)
			throw std::runtime_error("row or column isn't 3.");
	}
	catch (std::runtime_error err){
		return this->Determinant();
	}

	T add1 = element[0][0] * element[1][1] * element[2][2];
	T add2 = element[0][1] * element[1][2] * element[2][0];
	T add3 = element[0][2] * element[1][0] * element[2][1];
	T min1 = element[0][0] * element[1][2] * element[2][1];
	T min2 = element[0][1] * element[1][0] * element[2][2];
	T min3 = element[0][2] * element[1][1] * element[2][0];

	return add1 + add2 + add3 - min1 - min2 - min3;
}

template <typename T>
Matrix<T> Matrix<T>::Adjoint() const
{
	try{
		if (this->row != 3 || this->column != 3)
			throw std::runtime_error("row or column isn't 3.");
	}
	catch (std::runtime_error err){
		return *this;
	}

	// 余子式矩阵
	T a11 = element[1][1] * element[2][2] - element[1][2] * element[2][1];
	T a12 = element[1][2] * element[2][0] - element[1][0] * element[2][2];
	T a13 = element[1][0] * element[2][1] - element[1][1] * element[2][0];
	T a21 = element[0][2] * element[2][1] - element[0][1] * element[2][2];
	T a22 = element[0][0] * element[2][2] - element[0][2] * element[2][0];
	T a23 = element[0][1] * element[2][0] - element[0][0] * element[2][1];
	T a31 = element[0][1] * element[1][2] - element[0][2] * element[1][1];
	T a32 = element[0][2] * element[1][0] - element[0][0] * element[1][2];
	T a33 = element[0][0] * element[1][1] - element[0][1] * element[1][0];

	Matrix<T> result(3, 3, 0);
	result.element[0][0] = a11;
	result.element[0][1] = a12;
	result.element[0][2] = a13;
	result.element[1][0] = a21;
	result.element[1][1] = a22;
	result.element[1][2] = a23;
	result.element[2][0] = a31;
	result.element[2][1] = a32;
	result.element[2][2] = a33;

	return result.Transposition();
}

template <typename T> 
template <typename G>
Matrix<G> Matrix<T>::Inverse(G& det) const
{
	Matrix<G> errMat(3, 3, -1.0);
	try{
		if (this->row != 3 || this->column != 3)
			throw std::runtime_error("row or column isn't 3.");
	}
	catch (std::runtime_error err){
		return errMat;
	}
	try{
		if (this->Determinant() == 0)
			throw std::runtime_error("the determinant is 0.");
	}
	catch(std::runtime_error err){
		return errMat;
	}
	Matrix<T> adjMat = this->Adjoint();
	det = static_cast<G>(this->Determinant());

	Matrix<G> invMat(3, 3, 0.0);
	invMat = adjMat / det;
	return invMat;
}

template <typename T>
Matrix<T> SetDCM(float psi, float theta, float phi)
{
	if (psi > 180.0f)
		psi -= 360.0f;
	if (psi < -180.0f)
		psi += 360.0f;
	if (theta > 90.0f)
		theta -= 360.0f;
	if (theta < -90.0f)
		theta += 360.0f;
	if (phi > 180.0f)
		phi -= 360.0f;
	if (phi < -180.0f)
		phi += 360.0f;

	Matrix<T> DCM(3, 3, 0);
	float rad_psi = DEG2RAD(psi);
	float rad_theta = DEG2RAD(theta);
	float rad_phi = DEG2RAD(phi);

	DCM.element[0][0] = cosf(rad_psi) * cosf(rad_theta);
	DCM.element[0][1] = sinf(rad_psi) * cosf(rad_theta);
	DCM.element[0][2] = -1.0f * sinf(rad_theta);
	DCM.element[1][0] = cosf(rad_psi) * sinf(rad_theta) * sinf(rad_phi) - sinf(rad_psi) * cosf(rad_phi);
	DCM.element[1][1] = sinf(rad_psi) * sinf(rad_theta) * sinf(rad_phi) + cosf(rad_psi) * cosf(rad_phi);
	DCM.element[1][2] = cosf(rad_theta) * sinf(rad_phi);
	DCM.element[2][0] = cosf(rad_psi) * sinf(rad_theta) * cosf(rad_phi) + sinf(rad_psi) * sinf(rad_phi);
	DCM.element[2][1] = sinf(rad_psi) * sinf(rad_theta) * cosf(rad_phi) - cosf(rad_psi) * sinf(rad_phi);
	DCM.element[2][2] = cosf(rad_theta) * cosf(rad_phi);

	return DCM;
}

template <typename T>
Matrix<T> MakeMatrix(int row, int column, T num)
{
	Matrix<T> mat(row, column, num);
	return mat;
}

// 3D向量类
template <typename T>
class Vector3
{
public:
	Vector3();
	Vector3(T x, T y, T z);
	Vector3(const Vector3&);
	~Vector3();
	Vector3 operator=(const Vector3&);
	bool operator==(const Vector3&);
	bool operator!=(const Vector3&);
	void zero();								// 置为零向量
	Vector3 operator-(const Vector3&) const;
	Vector3 operator+(const Vector3&) const;
	template <typename G>
	Vector3<G> operator*(const G&) const;
	template <typename G>
	Vector3 operator/(const G&) const;
	Vector3 operator+=(const Vector3&);
	Vector3 operator-=(const Vector3&);
	template <typename G>
	Vector3 operator*=(const G&);
	template <typename G>
	Vector3 operator/=(const G&);
	void normalize();							// 向量标准化
	T operator*(const Vector3&) const;
	inline T mod()								// 求模
	{
		return sqrt(x * x + y * y + z * z);
	}
	Vector3 CrossMult(const Vector3&) const;	// 向量叉乘
	T x, y, z;
};

template <typename T>
Vector3<T>::Vector3()
{
	x = y = z = 0;
}

template <typename T>
Vector3<T>::Vector3(T x, T y, T z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}

template <typename T>
Vector3<T>::Vector3(const Vector3& vec)
{
	this->x = vec.x;
	this->y = vec.y;
	this->z = vec.z;
}

template <typename T>
Vector3<T>::~Vector3()
{
}

template <typename T>
Vector3<T> Vector3<T>::operator=(const Vector3& vec)
{
	x = vec.x;
	y = vec.y;
	z = vec.z;

	return *this;
}

template <typename T>
bool Vector3<T>::operator==(const Vector3& vec)
{
	return (x == vec.x && y == vec.y && z == vec.z);
}

template <typename T>
bool Vector3<T>::operator!=(const Vector3& vec)
{
	return (x != vec.x || y != vec.y || z != vec.z);
}

template <typename T>
void Vector3<T>::zero()
{
	x = y = z = 0;
}

template <typename T>
Vector3<T> Vector3<T>::operator-(const Vector3& vec) const
{
	Vector<T> result;
	result.x = x - vec.x;
	result.y = y - vec.y;
	result.z = z - vec.z;

	return result;
}

template <typename T>
Vector3<T> Vector3<T>::operator+(const Vector3& vec) const
{
	Vector<T> result;
	result.x = x + vec.x;
	result.y = y + vec.y;
	result.z = z + vec.z;

	return result;
}

template <typename T>
template <typename G>
Vector3<G> Vector3<T>::operator*(const G& scalar) const
{
	Vector3<G> result;
	result.x = x * scalar;
	result.y = y * scalar;
	result.z = z * scalar;

	return result;
}

template <typename T>
template <typename G>
Vector3<T> Vector3<T>::operator/(const G& scalar) const
{
	try{
		if (scalar == 0)
			throw std::runtime_error("scalar can't be zero.");
	}
	catch (std::runtime_error err){
		return *this;
	}

	Vector<T> result;
	result.x = x / scalar;
	result.y = y / scalar;
	result.z = z / scalar;

	return result;
}

template <typename T>
Vector3<T> Vector3<T>::operator+=(const Vector3& vec)
{
	x += vec.x;
	y += vec.y;
	z += vec.z;

	return *this;
}

template <typename T>
Vector3<T> Vector3<T>::operator-=(const Vector3& vec)
{
	x -= vec.x;
	y -= vec.y;
	z -= vec.z;

	return *this;
}

template <typename T>
template <typename G>
Vector3<T> Vector3<T>::operator*=(const G& scalar)
{
	x = x * scalar;
	y = y * scalar;
	z = z * scalar;
	
	return *this;
}

template <typename T>
template <typename G>
Vector3<T> Vector3<T>::operator/=(const G& scalar)
{
	x = x / scalar;
	y = y / scalar;
	z = z / scalar;
	
	return *this;
}

template <typename T>
void Vector3<T>::normalize()
{
	try{
		if (x == y == z == 0)
			throw std::runtime_error("vector can't be zero.");
	}
	catch (std::runtime_error err)
	{
		return *this;
	}

	T mod = sqrt(x * x + y * y + z * z);
	x /= mod;
	y /= mod;
	z /= mod;

	return *this;
}

template <typename T>
T Vector3<T>::operator*(const Vector3& vec) const
{
	return (x * vec.x + y * vec.y + z * vec.z);
}

template <typename T>
Vector3<T> Vector3 <T>::CrossMult(const Vector3& vec) const
{
	Vector3<T> result;
	result.x = y * vec.z - z * vec.y;
	result.y = z * vec.x - x * vec.z;
	result.z = x * vec.y - y * vec.x;

	return result;
}

template <typename T>
struct MyPoint
{
	T x;
	T y;
	T z;
};

// 判断多边形是凸多边形还是凹多边形
// 返回值: 0——凸多边形  1——凹多边形
template <typename T>
int PolygonIdentify(const MyPoint<T>* const array, const int& num)			// array: 顶点数组	num: 顶点个数
{
	Vector3<T> *vecarr = new Vector3<T>[num];
	for (int i = 0; i < num - 1; i++)
	{
		vecarr[i].x = array[i + 1].x - array[i].x;
		vecarr[i].y = array[i + 1].y - array[i].y;
		vecarr[i].z = 0;
	}
	vecarr[num - 1].x = array[0].x - array[num - 1].x;
	vecarr[num - 1].y = array[0].y - array[num - 1].y;
	vecarr[num - 1].z = 0;

	for (int i = 0; i < num - 1; i++)
	{
		if (vecarr[i] * vecarr[i + 1] < 0)
			return 1;
	}
	if (vecarr[num - 1] * vecarr[0] < 0)
		return 1;

	return 0;
}

// 将多边形分割成三角形
// 返回值: 顶点数组和顶点个数的pair
template <typename T>
std::pair<MyPoint<T>*, int> PolygonSlicer(const MyPoint<T>* const vertexarray, const int& num)
{
	std::pair<MyPoint<T>*, int> pairArray;
	int arrnum = 3 * (num - 2);
	MyPoint<T> *vexArray = new MyPoint<T>[arrnum];
	for (int i = 0; i < arrnum - 2; i++)
	{
		if (i % 3 == 0)
		{
			vexArray[i].x = vertexarray[0].x;
			vexArray[i].y = vertexarray[0].y;
			vexArray[i].z = vertexarray[0].z;
			vexArray[i + 1].x = vertexarray[i / 3 + 1].x;
			vexArray[i + 1].y = vertexarray[i / 3 + 1].y;
			vexArray[i + 1].z = vertexarray[i / 3 + 1].z;
			vexArray[i + 2].x = vertexarray[i / 3 + 2].x;
			vexArray[i + 2].y = vertexarray[i / 3 + 2].y;
			vexArray[i + 2].z = vertexarray[i / 3 + 2].z;
		}
	}
	pairArray.first = vexArray;
	pairArray.second = arrnum;

	return pairArray;
}

// quad q = [w, n] = [cos(a), n*sin(a)] = [cos(a) x*sin(a) y*sin(a) z*sin(a)]
template <typename T>
struct QUAD
{
    T w;
    T x;
    T y;
    T z;
}Quad;

Quad<T> quadLog(const Quad<T> *src)
{
    T alpha = acos(src->w);
    T multi = alpha / sinf(alpha);
    Quad<T> result;
    result.w = 0;
    result.x = src->x * multi;
    result.y = src->y * multi;
    result.z = src->z * multi;

    return result;
}

Quad<T> quadExp(const Quad<T> *src, float exponent)
{
    if(fabs(src->w) < 0.99999f)
    {
        T alpha = acos(src->w);
        T alpha_new = alpha * exponent;
        T multi = sinf(alpha_new) / sinf(alpha);
        Quad<T> result;
        result.w = cosf(alpha_new);
        result.x = src->x * multi;
        result.y = src->y * multi;
        result.z = src->z * multi;

        return result;
    }
    else
    {
        return 0;
    }
}

Quad<T> quadSlerp(const Quad<T> *quadOne, const Quad<T> *quadTwo, double t)
{
    Quad<T> result;
    double cosOmega = quadOne->w * quadTwo->w + quadOne->x * quadTwo->x + \
                      quadOne->y * quadTwo->y + quadOne->z * quadTwo->z;
    Quad<T> minusOne;
    minusOne.w = quadOne->w;
    minusOne.x = quadOne->x;
    minusOne.y = quadOne->y;
    minusOne.z = quadOne->z;
    if(cosOmega < 0.0)
    {
        minusOne.w = -1.0 * quadOne->w;
        minusOne.x = -1.0 * quadOne->x;
        minusOne.y = -1.0 * quadOne->y;
        minusOne.z = -1.0 * quadOne->z;
        cosOmega = -1.0 * cosOmega;
    }

    double k0, k1;
    if(cosOmega > 0.999999)
    {
        k0 = 1.0 - t;
        k1 = t;
    }
    else
    {
        double sinOmega = sqrt(1.0 - cosOmega * cosOmega);
        double omega = atan2(sinOmega, cosOmega);
        double oneOverSinOmega = 1.0 / sinOmega;
        k0 = sin((1.0 -t) * omega) * oneOverSinOmega;
        k1 = sin(t * omega) * oneOverSinOmega;
    }

    result.w = minusOne.w * k0 + quadTwo->w * k1;
    result.x = minusOne.x * k0 + quadTwo->x * k1;
    result.y = minusOne.y * k0 + quadTwo->y * k1;
    result.z = minusOne.z * k0 + quadTwo->z * k1;

    return result;
}
