#pragma once

#include <math.h>
#include <exception>
#include <stdexcept>
#include <vector>
#include <algorithm>

#define PI	3.1415926535897932384626
//#include "coordcnv.h"

//#ifndef DLL_API
//#define DECLDIR __declspec(dllimport)
//#else
//#define DECLDIR __declspec(dllexport)
//#endif

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
//DECLDIR double Angle(const double &destLat, const double &destLon, const POSITIONSTRUCT *P3DData);

// 线性变换算法
//DECLDIR double LinearTransform(const double& input, const double& InStart, const double& InEnd, const double& OutStart, const double& OutEnd);
namespace fan
{
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
    Matrix<T> Adjoint() const;					// 标准伴随矩阵/共轭矩阵
    template <typename G>
    Matrix<G> Inverse() const;					// 逆矩阵
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
            result.element[i][j] = this->element[i][j] + mat.element[i][j];
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
            result.element[i][j] = this->element[i][j] - mat.element[i][j];
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
Matrix<T> minorMatrix(Matrix<T> mat, int row, int column)
{
	Matrix<T> result(mat.row - 1, mat.column - 1, 0);
	int x = 0, y = 0;
	for (int i = 0; i < mat.row; i++)
	{
		for (int j = 0; j < mat.column; j++)
		{
			if (i != row && j != column)
			{
				result.element[x][y++] = mat.element[i][j];
				if (y == mat.column - 1)
				{
					x++;
					y = 0;
				}
			}
		}
	}
	return result;
}

template <typename T>
T Matrix<T>::Determinant() const
{
	try {
		if (this->row != this->column)
			throw std::runtime_error("this matrix is not a square matrix.");
	}
	catch (std::runtime_error) {
		return this->Determinant();
	}

	T result = 0;
	if (row == 1)
		result = this->element[0][0];
	else if (row == 2)
		result = this->element[0][0] * this->element[1][1] - this->element[0][1] * this->element[1][0];
	else
	{
		for (int i = 0; i < this->column; i++)
		{
			Matrix<T> temp = minorMatrix<T>(*this, 0, i);
			result += this->element[0][i] * pow(-1.0, i) * temp.Determinant();
		}
	}
	return result;
}

template <typename T>
Matrix<T> Matrix<T>::Adjoint() const
{
	try {
		if (this->row != this->column)
			throw std::runtime_error("this matrix is not a square matrix.");
	}
	catch (std::runtime_error) {
		return *this;
	}

	Matrix<T> result(this->row, this->column, 0);
	for (int i = 0; i < result.row; i++)
	{
		for (int j = 0; j < result.column; j++)
		{
			Matrix<T> temp = minorMatrix<T>(*this, i, j);
			result.element[i][j] = pow(-1.0, i + j) * temp.Determinant();
		}
	}
	return result.Transposition();
    // 余子式矩阵
    //T a11 = element[1][1] * element[2][2] - element[1][2] * element[2][1];
    //T a12 = element[1][2] * element[2][0] - element[1][0] * element[2][2];
    //T a13 = element[1][0] * element[2][1] - element[1][1] * element[2][0];
    //T a21 = element[0][2] * element[2][1] - element[0][1] * element[2][2];
    //T a22 = element[0][0] * element[2][2] - element[0][2] * element[2][0];
    //T a23 = element[0][1] * element[2][0] - element[0][0] * element[2][1];
    //T a31 = element[0][1] * element[1][2] - element[0][2] * element[1][1];
    //T a32 = element[0][2] * element[1][0] - element[0][0] * element[1][2];
    //T a33 = element[0][0] * element[1][1] - element[0][1] * element[1][0];

    //Matrix<T> result(3, 3, 0);
    //result.element[0][0] = a11;
    //result.element[0][1] = a12;
    //result.element[0][2] = a13;
    //result.element[1][0] = a21;
    //result.element[1][1] = a22;
    //result.element[1][2] = a23;
    //result.element[2][0] = a31;
    //result.element[2][1] = a32;
    //result.element[2][2] = a33;

    //return result.Transposition();
}

template <typename T>
template <typename G>
Matrix<G> Matrix<T>::Inverse() const
{
	try {
		if (this->row != this->column)
			throw std::runtime_error("this matrix is not a square matrix.");
	}
	catch (std::runtime_error) {
		return Matrix<G>(1, 1, 0);
	}

	Matrix<T> adj = this->Adjoint();
	Matrix<G> result(this->row, this->column, 0);
	T det = this->Determinant();
	for (int i = 0; i < result.row; i++)
	{
		for (int j = 0; j < result.column; j++)
		{
			result.element[i][j] = (G)adj.element[i][j] / det;
		}
	}
	return result;
}

template <typename T>
Matrix<T> SetDCM(float psi, float theta, float phi)		// h r p
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
    float rad_psi = DEG2RAD(-psi);
    float rad_theta = DEG2RAD(-theta);
    float rad_phi = DEG2RAD(-phi);

	//DCM.element[0][0] = cosf(rad_psi) * cosf(rad_theta);
	//DCM.element[0][1] = sinf(rad_phi) * sinf(rad_theta) * cosf(rad_psi) - sinf(rad_psi) * cosf(rad_phi);
	//DCM.element[0][2] = cosf(rad_psi) * sinf(rad_theta) * cosf(rad_phi) + sinf(rad_psi) * sinf(rad_phi);
	//DCM.element[1][0] = sinf(rad_phi) * cosf(rad_theta);
	////DCM.element[1][0] = sinf(rad_psi) * cosf(rad_theta);
	//DCM.element[1][1] = sinf(rad_psi) * sinf(rad_theta) * sinf(rad_phi) + cosf(rad_psi) * cosf(rad_phi);
	//DCM.element[1][2] = sinf(rad_psi) * sinf(rad_theta) * cosf(rad_phi) - cosf(rad_psi) * sinf(rad_phi);
	//DCM.element[2][0] = -1.0f * sinf(rad_theta);
	//DCM.element[2][1] = cosf(rad_theta) * sinf(rad_phi);
	//DCM.element[2][2] = cosf(rad_theta) * cosf(rad_phi);
	Matrix<T> RotX(3, 3, 0);
	Matrix<T> RotY(3, 3, 0);
	Matrix<T> RotZ(3, 3, 0);
	RotX.element[0][0] = 1.0;
	RotX.element[0][1] = 0.0;
	RotX.element[0][2] = 0.0;
	RotX.element[1][0] = 0.0;
	RotX.element[1][1] = cos(rad_theta);
	RotX.element[1][2] = sin(rad_theta);
	RotX.element[2][0] = 0.0;
	RotX.element[2][1] = -1.0f * sin(rad_theta);
	RotX.element[2][2] = cos(rad_theta);

	RotY.element[0][0] = cos(rad_phi);
	RotY.element[0][1] = 0.0;
	RotY.element[0][2] = -1.0f * sin(rad_phi);
	RotY.element[1][0] = 0.0;
	RotY.element[1][1] = 1.0;
	RotY.element[1][2] = 0.0;
	RotY.element[2][0] = sin(rad_phi);
	RotY.element[2][1] = 0;
	RotY.element[2][2] = cos(rad_phi);

	RotZ.element[0][0] = cos(rad_psi);
	RotZ.element[0][1] = sin(rad_psi);
	RotZ.element[0][2] = 0.0;
	RotZ.element[1][0] = -1.0f * sin(rad_psi);
	RotZ.element[1][1] = cos(rad_psi);
	RotZ.element[1][2] = 0.0;
	RotZ.element[2][0] = 0.0;
	RotZ.element[2][1] = 0;
	RotZ.element[2][2] = 1.0;
	// Z->X->Y
	Matrix<T> temp = RotZ * RotX;
	DCM = temp * RotY;
    return DCM;
}

template <typename T>
Matrix<T> MakeMatrix(int row, int column, T num)
{
    Matrix<T> mat(row, column, num);
    return mat;
}

template <typename T>
Matrix<T> matrixMul(Matrix<T> one, Matrix<T> two)
{
	Matrix<T> result(one.row, two.column, 0);
    for (int i = 0; i < one.row; i++)
    {
        for (int j = 0; j < two.column; j++)
        {
            for (int k = 0; k < one.column; k++)
                result.element[i][j] += one.element[i][k] * two.element[k][j];
        }
    }
    return result;
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
    Vector3 normalize();							// 向量标准化
    T operator*(const Vector3<T>&) const;
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
    Vector3<T> result;
    result.x = x - vec.x;
    result.y = y - vec.y;
    result.z = z - vec.z;

    return result;
}

template <typename T>
Vector3<T> Vector3<T>::operator+(const Vector3& vec) const
{
    Vector3<T> result;
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

    Vector3<T> result;
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
Vector3<T> Vector3<T>::normalize()
{
    try{
        if (x == 0 && y ==0 && z == 0)
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
T Vector3<T>::operator*(const Vector3<T>& vec) const
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

// 2D向量类
template <typename T>
class Vector2
{
public:
	Vector2();
	Vector2(T x, T y);
	Vector2(const Vector2&);
	~Vector2();
	Vector2 operator=(const Vector2&);
	bool operator==(const Vector2&);
	bool operator!=(const Vector2&);
	void zero();								// 置为零向量
	Vector2 operator-(const Vector2&) const;
	Vector2 operator+(const Vector2&) const;
	template <typename G>
	Vector2<G> operator*(const G&);
	template <typename G>
	Vector2 operator/(const G&) const;
	Vector2 operator+=(const Vector2&);
	Vector2 operator-=(const Vector2&);
	template <typename G>
	Vector2 operator*=(const G&);
	template <typename G>
	Vector2 operator/=(const G&);
	Vector2 normalize();							// 向量标准化
	T operator*(const Vector2<T>&);
	inline T mod()								// 求模
	{
		return sqrt(x * x + y * y);
	}
	Vector2 CrossMult(const Vector2&) const;	// 向量叉乘
	T x, y;
};

template<typename T>
Vector2<T>::Vector2()
{
	x = y = 0;
}

template<typename T>
Vector2<T>::Vector2(T x, T y)
{
	this->x = x;
	this->y = y;
}

template <typename T>
Vector2<T>::Vector2(const Vector2& vec)
{
	this->x = vec.x;
	this->y = vec.y;
}

template <typename T>
Vector2<T>::~Vector2()
{

}

template<typename T>
Vector2<T> Vector2<T>::operator =(const Vector2<T>& vec)
{
	x = vec.x;
	y = vec.y;
	return *this;
}

template<typename T>
T Vector2<T>::operator *(const Vector2<T>& vec)
{
	return (x * vec.x + y * vec.y);
}

template <typename T>
template <typename G>
Vector2<G> Vector2<T>::operator *(const G& Scalar)
{
	x = x * Scalar;
	y = y * Scalar;
	return *this;
}

template<typename T>
template<typename G>
Vector2<T> Vector2<T>::operator *=(const G& Scalar)
{
	x *= Scalar;
	y *= Scalar;
	return *this;
}

template<typename T>
Vector2<T> Vector2<T>::normalize()
{
	if(x != 0 && y != 0)
	{
		x = x / mod();
		y = y / mod();
	}

	return *this;
}


// vector to Euler angle. roll = 0.0
template <typename T>
void vec3_To_Euler(const Vector3<T>& vec3, T Euler[3])
{
	if(vec3.y > 0.0)
	{
		Euler[0] = atan(vec3.x / vec3.y);
		if(vec3.x == 0.0)
			Euler[0] = 0.0;
	}
	if (vec3.y < 0.0)
	{
		Euler[0] = -1.0 * PI + abs(vec3.x) / vec3.x * atan(vec3.x / vec3.y);
	}
	if(vec3.y == 0.0)
	{
		Euler[0] = -1.0 * PI / 2.0 * abs(vec3.x) / vec3.x;
	}
	if(vec3.x !=0 || vec3.y != 0)
		Euler[1] = atan(vec3.z / sqrt(vec3.x * vec3.x + vec3.y * vec3.y));
	else
		Euler[1] = PI * abs(vec3.z) / vec3.z;
	Euler[0] = RAD2DEG(Euler[0]);
	Euler[1] = RAD2DEG(Euler[1]);
	Euler[2] = 0.0;
}

template <typename T>
struct MyPoint
{
    T x;
    T y;
    T z;
    MyPoint<T> operator =(const MyPoint<T>&);
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

// QUAd q = [w, n] = [cos(a), n*sin(a)] = [cos(a) x*sin(a) y*sin(a) z*sin(a)]
template <typename T>
struct QUAd
{
    T w;
    T x;
    T y;
    T z;
};

template <typename T>
QUAd<T> QUAdLog(const QUAd<T> *src)
{
    T alpha = acos(src->w);
    T multi = alpha / sinf(alpha);
    QUAd<T> result;
    result.w = 0;
    result.x = src->x * multi;
    result.y = src->y * multi;
    result.z = src->z * multi;

    return result;
}

template <typename T>
QUAd<T> QUAdExp(const QUAd<T> *src, float exponent)
{
    if(fabs(src->w) < 0.99999f)
    {
        T alpha = acos(src->w);
        T alpha_new = alpha * exponent;
        T multi = sinf(alpha_new) / sinf(alpha);
        QUAd<T> result;
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

template <typename T>
QUAd<T> QUAdSlerp(const QUAd<T> *QUAdOne, const QUAd<T> *QUAdTwo, double t)
{
    QUAd<T> result;
    double cosOmega = QUAdOne->w * QUAdTwo->w + QUAdOne->x * QUAdTwo->x + \
                      QUAdOne->y * QUAdTwo->y + QUAdOne->z * QUAdTwo->z;
    QUAd<T> minusOne;
    minusOne.w = QUAdOne->w;
    minusOne.x = QUAdOne->x;
    minusOne.y = QUAdOne->y;
    minusOne.z = QUAdOne->z;
    if(cosOmega < 0.0)
    {
        minusOne.w = -1.0 * QUAdOne->w;
        minusOne.x = -1.0 * QUAdOne->x;
        minusOne.y = -1.0 * QUAdOne->y;
        minusOne.z = -1.0 * QUAdOne->z;
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

    result.w = minusOne.w * k0 + QUAdTwo->w * k1;
    result.x = minusOne.x * k0 + QUAdTwo->x * k1;
    result.y = minusOne.y * k0 + QUAdTwo->y * k1;
    result.z = minusOne.z * k0 + QUAdTwo->z * k1;

    return result;
}

/*
说明:
    MyPoint: 点结构体
    sortX: 将点以X坐标判断大小
    sortY: 将点以Y坐标判断大小
    LineSegment: 线段结构体
    makeLineSegment: 两点组成线段
    lineSegmentLength: 求线段长
    IsCross: 判断两线段是否有交点
    QUAdrangle: 四边形结构体
    makeQUAdrangle: 将四点组成四边形
    IsConvexQUAdr: 判断是否标准四边形
    makeStandardQUAd: 将四边形顶点按标准方式排序
    halfLingSegment: 将线段两等分
    quarterLineSegment: 将线段四等分
        PtInQUAd: 判断点是否在矩形内
*/


//template <typename T>
//struct MyPoint
//{
//	T x;
//	T y;
//	MyPoint<T> operator=(const MyPoint<T>& source);
//};

template <typename T>
MyPoint<T> MyPoint<T>::operator=(const MyPoint<T>& source)
{
    this->x = source.x;
    this->y = source.y;
    this->z = source.z;
    return *this;
}

template <typename T>
bool sortX(MyPoint<T> pt1, MyPoint<T> pt2)
{
    return pt1.x < pt2.x;
}

template <typename T>
bool sortY(MyPoint<T> pt1, MyPoint<T> pt2)
{
    return pt1.y < pt2.y;
}

template <typename T>
struct LineSegment
{
    MyPoint<T> P1;
    MyPoint<T> P2;
};

template <typename T>
LineSegment<T> makeLineSegment(MyPoint<T> A, MyPoint<T> B)
{
    LineSegment<T> result;
    result.P1.x = A.x;
    result.P1.y = A.y;
    result.P2.x = B.x;
    result.P2.y = B.y;

    return result;
}

template <typename T>
double lineSegmentLength(LineSegment<T> L1)
{
    return sqrt((L1.P1.x - L1.P2.x) * (L1.P1.x - L1.P2.x) + (L1.P1.y - L1.P2.y) * (L1.P1.y - L1.P2.y));
}

template <typename T>
bool IsCross(LineSegment<T> AB, LineSegment<T> CD)
{
    double kAB, bAB, kCD, bCD;

    if (AB.P2.x != AB.P1.x)
        kAB = (double)(AB.P2.y - AB.P1.y) / (double)(AB.P2.x - AB.P1.x);
    else
        return false;
    bAB = AB.P1.y - kAB * AB.P1.x;
    if (CD.P2.x != CD.P1.x)
        kCD = (double)(CD.P2.y - CD.P1.y) / (double)(CD.P2.x - CD.P1.x);
    else
        return false;

    if (kAB == kCD)
    {
        bCD = CD.P1.y - kCD * CD.P1.x;
        if (bAB == bCD)
            return false;
        else
            return true;
    }

    LineSegment<T> AC = makeLineSegment(AB.P1, CD.P1);
    LineSegment<T> BD = makeLineSegment(AB.P2, CD.P2);

    T minACX = AC.P1.x < AC.P2.x ? AC.P1.x : AC.P2.x;
    T maxACX = AC.P1.x > AC.P2.x ? AC.P1.x : AC.P2.x;
    T minACY = AC.P1.y < AC.P2.y ? AC.P1.y : AC.P2.y;
    T maxACY = AC.P1.y > AC.P2.y ? AC.P1.y : AC.P2.y;
    T minBDX = BD.P1.x < BD.P2.x ? BD.P1.x : BD.P2.x;
    T maxBDX = BD.P1.x > BD.P2.x ? BD.P1.x : BD.P2.x;
    T minBDY = BD.P1.y < BD.P2.y ? BD.P1.y : BD.P2.y;
    T maxBDY = BD.P1.y > BD.P2.y ? BD.P1.y : BD.P2.y;

    double kAC, bAC, kBD, bBD, crossX, crossY;
    if (AC.P2.x != AC.P1.x)
    {
        kAC = (double)(AC.P2.y - AC.P1.y) / (double)(AC.P2.x - AC.P1.x);
        bAC = (double)(AC.P1.y - kAC * AC.P1.x);
    }
    else
    {
        kBD = (double)(BD.P2.y - BD.P1.y) / (double)(BD.P2.x - BD.P1.x);
        bBD = BD.P1.y - kBD * BD.P1.x;

        crossX = AC.P1.x;
        crossY = kBD * crossX + bBD;
    }

    if (BD.P2.x != BD.P1.x)
    {
        kBD = (double)(BD.P2.y - BD.P1.y) / (double)(BD.P2.x - BD.P1.x);
        bBD = BD.P1.y - kBD * BD.P1.x;
    }
    else
    {
        kAC = (double)(AC.P2.y - AC.P1.y) / (double)(AC.P2.x - AC.P1.x);
        bAC = (AC.P1.y - kAC * AC.P1.x);
        crossX = BD.P1.x;
        crossY = kAC * crossX + bAC;
    }

    if (AC.P2.x != AC.P1.x && BD.P2.x != BD.P1.x)
    {
        crossX = (double)(bBD - bAC) / (double)(kAC - kBD);
        crossY = kAC * crossX + bAC;
    }

    if (crossX >= minACX && crossX <= maxACX && crossX >= minBDX && crossX <= maxBDX && \
        crossY >= minACY && crossY <= maxACY && crossY >= minBDY && crossY <= maxBDY)
        return true;
    else
        return false;
}

template <typename T>
struct QUAdrangle
{
    MyPoint<T> A;
    MyPoint<T> B;
    MyPoint<T> C;
    MyPoint<T> D;
};

template <typename T>
QUAdrangle<T> makeQUAdrangle(std::vector<MyPoint<T>> vecQUAd)
{
    QUAdrangle<T> result;
    result.A.x = vecQUAd[0].x;
    result.A.y = vecQUAd[0].y;
    result.B.x = vecQUAd[1].x;
    result.B.y = vecQUAd[1].y;
    result.C.x = vecQUAd[2].x;
    result.C.y = vecQUAd[2].y;
    result.D.x = vecQUAd[3].x;
    result.D.y = vecQUAd[3].y;

    return result;
}

template <typename T>
bool IsConvexQUAdr(QUAdrangle<T> QUAd)
{
    LineSegment<T> AB = makeLineSegment(QUAd.A, QUAd.B);
    LineSegment<T> CD = makeLineSegment(QUAd.C, QUAd.D);

    return IsCross(AB, CD);
}

template <typename T>
QUAdrangle<T> makeStandardQUAd(QUAdrangle<T> source)
{
//    QUAdrangle<T> result;
//	std::vector<MyPoint<T>> vecPt;
//	vecPt.push_back(source.A);
//	vecPt.push_back(source.B);
//	vecPt.push_back(source.C);
//	vecPt.push_back(source.D);

//	std::sort(vecPt.begin(), vecPt.end(), sortX<T>);
//	std::vector<MyPoint<T>>::iterator iterVec;
//	iterVec = vecPt.begin();
//	result.A = iterVec->y < (iterVec + 1)->y ? *iterVec : *(iterVec + 1);
//	result.B = iterVec->y > (iterVec + 1)->y ? *iterVec : *(iterVec + 1);
//	result.C = (iterVec + 2)->y > (iterVec + 3)->y ? *(iterVec + 2) : *(iterVec + 3);
//	result.D = (iterVec + 2)->y < (iterVec + 3)->y ? *(iterVec + 2) : *(iterVec + 3);

//	return result;
}

template <typename G, typename T>
std::vector<MyPoint<G>> halfLingSegment(LineSegment<T> source)
{
    MyPoint<G> halfPoint;
    std::vector<MyPoint<G>> result;
    MyPoint<G> _vector;

    _vector.x = static_cast<G>(source.P1.x);
    _vector.y = static_cast<G>(source.P1.y);

    result.push_back(_vector);

    _vector.x = source.P2.x - source.P1.x;
    _vector.y = source.P2.y - source.P1.y;
    _vector.x /= 2.0;
    _vector.y /= 2.0;

    halfPoint.x = source.P1.x + _vector.x;
    halfPoint.y = source.P1.y + _vector.y;

    result.push_back(halfPoint);
    _vector.x = static_cast<G>(source.P2.x);
    _vector.y = static_cast<G>(source.P2.y);
    result.push_back(_vector);

    return result;
}


template <typename G, typename T>
std::vector<MyPoint<G>> quarterLineSegment(LineSegment<T> source)
{
    MyPoint<G> pushPoint;
    std::vector<MyPoint<G>> result;
    MyPoint<G> _vector;

    _vector.x = static_cast<G>(source.P1.x);
    _vector.y = static_cast<G>(source.P1.y);
    result.push_back(_vector);

    _vector.x = source.P2.x - source.P1.x;
    _vector.y = source.P2.y - source.P1.y;

    MyPoint<double> temp = _vector;

    for (double para = 0.25; para < 1.0; para += 0.25)
    {
        temp.x = _vector.x * para;
        temp.y = _vector.y * para;
        pushPoint.x = source.P1.x + temp.x;
        pushPoint.y = source.P1.y + temp.y;
        result.push_back(pushPoint);
    }

    _vector.x = static_cast<G>(source.P2.x);
    _vector.y = static_cast<G>(source.P2.y);
    result.push_back(_vector);

    return result;
}

template <typename T, typename G>
bool PtInQUAd(MyPoint<T> pt, QUAdrangle<G> QUAd)
{
    MyPoint<G> _vecAB;
    _vecAB.x = QUAd.B.x - QUAd.A.x;
    _vecAB.y = QUAd.B.y - QUAd.A.y;

    MyPoint<G> _vecBC;
    _vecBC.x = QUAd.C.x - QUAd.B.x;
    _vecBC.y = QUAd.C.y - QUAd.B.y;

    MyPoint<G> _vecCD;
    _vecCD.x = QUAd.D.x - QUAd.C.x;
    _vecCD.y = QUAd.D.y - QUAd.C.y;

    MyPoint<G> _vecDA;
    _vecDA.x = QUAd.A.x - QUAd.D.x;
    _vecDA.y = QUAd.A.y - QUAd.D.y;

    MyPoint<double> ABP;
    ABP.y = pt.y;
    ABP.x = QUAd.A.x + _vecAB.x * (pt.y - QUAd.A.y) / _vecAB.y;

    MyPoint<double> BCP;
    BCP.x = pt.x;
    BCP.y = _vecBC.y / _vecBC.x * (pt.x - QUAd.B.x) + QUAd.B.y;

    MyPoint<double> CDP;
    CDP.y = pt.y;
    CDP.x = _vecCD.x / _vecCD.y * (pt.y - QUAd.C.y) + QUAd.C.x;

    MyPoint<double> DAP;
    DAP.x = pt.x;
    DAP.y = _vecDA.y / _vecDA.x * (pt.x - QUAd.D.x) + QUAd.D.y;

    return (pt.x > ABP.x) && (pt.y < BCP.y) && (pt.x < CDP.x) && (pt.y > DAP.y);
}

template <typename T>
Matrix<T> RodrigueMatrix(Vector3<T> src, Vector3<T> dst)
{
    double arctheta = (src * dst) / src.mod() / dst.mod();
    double theta = acosf(arctheta);
    Matrix<T> result(3, 3, 0);
    result.element[0][0] = cos(theta) + src.x * src.x * (1 - cos(theta));
    result.element[0][1] = src.x * src.y * (1 - cos(theta)) - src.z * sin(theta);
    result.element[0][2] = src.y * sin(theta) + src.x * src.z * (1 - cos(theta));
    result.element[1][0] = src.z * sin(theta) + src.x * src.y * (1 - cos(theta));
    result.element[1][1] = cos(theta) + src.y * src.y * (1 - cos(theta));
    result.element[1][2] = -1.0 * src.x * sin(theta) + src.y * src.z * (1 - cos(theta));
    result.element[2][0] = -1.0 * src.y * sin(theta) + src.x * src.z * (1 - cos(theta));
    result.element[2][1] = src.x * sin(theta) + src.y * src.z * (1 - cos(theta));
    result.element[2][2] = cos(theta) + src.z * src.z * (1 - cos(theta));
}

// myTrackTrail: 两点追踪渐近
// DstPos[3]: 目标点位置
// DstRot[3]: 目标点朝向(暂时无用)
// SrcPos[3]: 前进点位置
// SrcRot[3]: 前进点朝向
// CloseEnough: 判定有效距离
// Speed:		前进点行进速度

template<typename T>
bool myTrackTrail(T DstPos[3], T DstRot[3], T SrcPos[3], T SrcRot[3], T CloseEnough, T Speed)
{
	T deltaX = DstPos[0] - SrcPos[0];
	T deltaY = DstPos[1] - SrcPos[1];
	T deltaZ = DstPos[2] - SrcPos[2];
	//double deltaH = RAD2DEG(atan(deltaX / deltaY))/* - m_missilePos.yaw*/;
	T deltaH;
	if(abs(deltaY - 0.0) < 0.00001)
	{
		if(deltaX > 0.0)
		{
			deltaH = -90.0;
		}
		if(deltaX < 0.0)
			deltaH = 90.0;
		if(deltaX - 0.0 < 0.00001)
			deltaH = 0.0;
	}
	if(abs(deltaX - 0.0) < 0.0001)
	{
		if(deltaY > 0.0)
			deltaH = 0.0;
		if(deltaY < 0.0)
			deltaH = 180.0;
	}
	if(deltaY > 0.0 && !(abs(deltaX - 0.0) < 0.0001))
	{
		deltaH = -1.0 * RAD2DEG(atan(deltaX / deltaY));
	}
	if(deltaY < 0.0 && deltaX < 0.0)
	{
		deltaH = 180.0 - RAD2DEG(atan(abs(deltaX) / abs(deltaY)));
	}
	if(deltaY < 0.0 && deltaX > 0.0)
	{
		deltaH = -180.0 + RAD2DEG(atan(deltaX / abs(deltaY)));
	}
	double deltaP = RAD2DEG(asin(deltaZ / sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ)))/* - m_missilePos.pitch*/;

	SrcRot[0] = deltaH;
	SrcRot[1] = deltaP;
	SrcRot[2] = 0.0;

	Matrix<T> m_DCM = SetDCM<T>(-SrcRot[0], -SrcRot[1], -SrcRot[2]);
	m_DCM.Transposition();
	Matrix<T> m_pos(1, 3, 0);
	m_pos.element[0][0] = 0.0;
	m_pos.element[0][1] = Speed;
	m_pos.element[0][2] = 0.0;
	Matrix<T> dst(1, 3, 0);
	dst = matrixMul<T>(m_pos, m_DCM);
	T tranX, tranY, tranZ;
	tranX = dst.element[0][0];
	tranY = dst.element[0][1];
	tranZ = dst.element[0][2];
	SrcPos[0] += tranX;
	SrcPos[1] += tranY;
	SrcPos[2] += tranZ;

	if(deltaX *deltaX + deltaY * deltaY + deltaZ * deltaZ < CloseEnough * CloseEnough)
	{
		return true;
	}
	return false;
}

template <typename T>
bool myTrackTrail2(T DstPos[3], T DstRot[3], T SrcPos[3], T SrcRot[3], T CloseEnough, T Speed)
{
	Vector2<T> vecA(sin(DEG2RAD(DstRot[0])), cos(DEG2RAD(DstRot[0])));
	Vector2<T> vecAB(SrcPos[0]- DstPos[0], SrcPos[1] - DstPos[1]);
	T cosTheta = vecA * vecAB / (vecA.mod() * vecAB.mod());
	vecA.normalize();
	T vecTemp = vecAB.mod() * cosTheta;
	Vector2<T> vecAD = vecA.operator *<T>(vecTemp);
	vecAD *= 0.65;
	double D[3];
	D[0] = DstPos[0] + vecAD.x;
	D[1] = DstPos[1] + vecAD.y;
	D[2] = DstPos[2] + 0.0;

	return myTrackTrail<double>(D, NULL, SrcPos, SrcRot, CloseEnough, Speed);
}

template <typename T>
void worldToScreen(T eye[3], T euler[3], T dst[3], T viewport[2], T* xy)
{
	Matrix<T> DCM = SetDCM<T>(euler[0], euler[2], euler[1]);
	DCM.Transposition();
	Matrix<T>  mat(1, 3, 0);
	mat.element[0][1] = 1.0;
	Matrix<T> Odot(1, 3, 0);
	Odot = matrixMul<T>(mat, DCM);
	Vector3<T> OOdot(Odot.element[0][0], Odot.element[0][1], Odot.element[0][2]);
	Vector3<T> OA(dst[0] - eye[0], dst[1] - eye[1], dst[2] - eye[2]);
	T cosAlpha = OOdot * OA / OA.mod() / OOdot.mod();
	T modOD = OA.mod() * cosAlpha;
	Vector3<T> OD = OOdot * modOD;
	Vector3<T> M(OD.x + eye[0], OD.y + eye[1], OD.z + eye[2]);
	T AB = 2 * modOD * tan(DEG2RAD(viewport[1] / 2.0));
	T BC = 2 * modOD * tan(DEG2RAD(viewport[0] / 2.0));
	//T down[3] = {euler[0], euler[1] - 90.0, euler[2]};
	//mat.element[0][1] = AB / 2;
	//DCM = SetDCM<T>(down[0], down[2], down[1]);
	//DCM.Transposition();
	//Matrix<T> ME(1, 3, 0);
	//ME = matrixMul<T>(mat, DCM);
	//Vector3<T> E(ME.element[0][0] + M.x, ME.element[0][1] + M.y, ME.element[0][2] + M.z);
	//T left[3] = { euler[0] + 90.0, euler[1], euler[2] };
	//mat.element[0][1] = BC / 2;
	//DCM = SetDCM<T>(left[0], left[2], left[1]);
	//DCM.Transposition();
	//Matrix<T> EB(1, 3, 0);
	//EB = matrixMul<T>(mat, DCM);
	//Vector3<T> B(EB.element[0][0] + E.x, EB.element[0][1] + E.y, EB.element[0][2] + E.z);
	//T right[3] = { euler[0] - 90.0, euler[1], euler[2] };
	//mat.element[0][1] = BC;
	//DCM = SetDCM<T>(right[0], right[2], right[1]);
	//DCM.Transposition();
	//Matrix<T> CB(1, 3, 0);
	//CB = matrixMul<T>(mat, DCM);
	//Vector3<T> C(CB.element[0][0] + B.x, CB.element[0][1] + B.y, CB.element[0][2] + B.z);
	//Vector3<T> vecBC(C.x - B.x, C.y - B.y, C.z - B.z);
	//T up[3] = { euler[0], euler[1] + 90.0, euler[2] };
	//mat.element[0][1] = AB;
	//DCM = SetDCM<T>(up[0], up[2], up[1]);
	//DCM.Transposition();
	//Matrix<T> BA(1, 3, 0);
	//BA = matrixMul<T>(mat, DCM);
	//Vector3<T> A(BA.element[0][0] + B.x, BA.element[0][1] + B.y, BA.element[0][2] + B.z);
	//mat.element[0][1] = sqrt(OD.mod() * OD.mod() + 0.25f * AB * AB + 0.25f * BC * BC);
	//T theta = RAD2DEG(asin(AB * 0.5 / mat.element[0][1]));
	//DCM = SetDCM<T>(euler[0] + 0.5f * viewport[0], euler[2], euler[1] - theta);
	//DCM.Transposition();
	//Matrix<T> matOB = matrixMul<T>(mat, DCM);
	//Vector3<T> B(matOB.element[0][0] + eye[0], matOB.element[0][1] + eye[1], matOB.element[0][2] + eye[2]);
	//DCM = SetDCM<T>(euler[0] + 0.5 * viewport[0], euler[2], euler[1] + theta);
	//DCM.Transposition();
	////mat.element[0][1] = sqrt(OD.mod() * OD.mod() + 0.25 * AB * AB + 0.25 * BC * BC);
	//Matrix<T> matOA = matrixMul<T>(mat, DCM);
	//Vector3<T> A(matOA.element[0][0] + eye[0], matOA.element[0][1] + eye[1], matOA.element[0][2] + eye[2]);
	//DCM = SetDCM<T>(euler[0] - 0.5 * viewport[0], euler[2], euler[1] - theta);
	//DCM.Transposition();
	////mat.element[0][1] = sqrt(OD.mod() * OD.mod() + 0.25 * AB * AB + 0.25 * BC * BC);
	//Matrix<T> matOC = matrixMul<T>(mat, DCM);
	//Vector3<T> C(matOC.element[0][0] + eye[0], matOC.element[0][1] + eye[1], matOC.element[0][2] + eye[2]);
	DCM = SetDCM<T>(euler[0], euler[2], euler[1] - viewport[1] * 0.5);
	DCM.Transposition();
	mat.element[0][1] = sqrt(OD.mod() * OD.mod() + 0.25f * AB * AB);
	Matrix<T> OE = matrixMul(mat, DCM);
	Vector3<T> E(OE.element[0][0] + eye[0], OE.element[0][1] + eye[1], OE.element[0][2] + eye[2]);
	Vector3<T> ME(E.x - M.x, E.y - M.y, E.z - M.z);

	DCM = SetDCM<T>(euler[0], euler[2], euler[1] + viewport[1] * 0.5);
	DCM.Transposition();
	mat.element[0][1] = sqrt(OD.mod() * OD.mod() + 0.25f * AB * AB);
	Matrix<T> OG = matrixMul(mat, DCM);
	Vector3<T> G(OG.element[0][0] + eye[0], OG.element[0][1] + eye[1], OG.element[0][2] + eye[2]);
	Vector3<T> MG(G.x - M.x, G.y - M.y, G.z - M.z);

	DCM = SetDCM<T>(euler[0] + viewport[0] * 0.5, euler[2], euler[1]);
	DCM.Transposition();
	mat.element[0][1] = sqrt(OD.mod() * OD.mod() + 0.25f * BC * BC);
	Matrix<T> OF = matrixMul(mat, DCM);
	Vector3<T> F(OF.element[0][0] + eye[0], OF.element[0][1] + eye[1], OF.element[0][2] + eye[2]);
	Vector3<T> MF(F.x - M.x, F.y - M.y, F.z - M.z);

	DCM = SetDCM<T>(euler[0] - viewport[0] * 0.5, euler[2], euler[1]);
	DCM.Transposition();
	mat.element[0][1] = sqrt(OD.mod() * OD.mod() + 0.25f * BC * BC);
	Matrix<T> OH = matrixMul(mat, DCM);
	Vector3<T> H(OH.element[0][0] + eye[0], OH.element[0][1] + eye[1], OH.element[0][2] + eye[2]);
	Vector3<T> MH(H.x - M.x, H.y - M.y, H.z - M.z);

	Vector3<T> MB = ME + MF;
	Vector3<T> MA = MF + MG;
	Vector3<T> MC = ME + MH;
	Vector3<T> B(MB.x + M.x, MB.y + M.y, MB.z + M.z);
	Vector3<T> A(MA.x + M.x, MA.y + M.y, MA.z + M.z);
	Vector3<T> C(MC.x + M.x, MC.y + M.y, MC.z + M.z);

	Vector3<T> vecBC(C.x - B.x, C.y - B.y, C.z - B.z);
	Vector3<T> vecBA(A.x - B.x, A.y - B.y, A.z - B.z);
	Vector3<T> BN(dst[0] - B.x, dst[1] - B.y, dst[2] - B.z);

	T cosBeta = BN * vecBC / BN.mod() / vecBC.mod();
	xy[0] = BN.mod() * cosBeta / vecBC.mod();
	xy[1] = BN.mod() * sqrt(1 - cosBeta * cosBeta) / vecBC.mod();

	xy[0] = 2 * xy[0] - 1;
	xy[1] = 2 * xy[1] - 1;

	if (abs(xy[0]) < 0.000001)
		xy[0] = 0.0;
	if (abs(xy[1]) < 0.000001)
		xy[1] = 0.0;
}

template <typename T>
void worldToScreen_Matrix(T eyepos[3], T eyeEuler[3], T viewport[2], T objpos[3], T screen[2])
{
	Matrix<T> dcm = SetDCM<T>(eyeEuler[0], eyeEuler[1], eyeEuler[2]);
	dcm.Transposition();
	Matrix<T> conv(4, 4, 0);
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0;j < 3; j++)
		{
			conv.element[i][j] = dcm.element[i][j];
		}
	}
	conv.element[3][0] = eyepos[0];
	conv.element[3][1] = eyepos[1];
	conv.element[3][2] = eyepos[2];
	conv.element[3][3] = 1.0;
	conv = conv.Inverse<T>();
	Matrix<T> obj(1, 4, 1);
	for (int i = 0; i < 3; i++)
	{
		obj.element[0][i] = objpos[i];
	}
	obj = obj * conv;
	T LimitX = obj.element[0][1] * tanf(DEG2RAD(viewport[0] / 2.0f));
	T LimitY = obj.element[0][1] * tanf(DEG2RAD(viewport[1] / 2.0f));
	screen[0] = obj.element[0][0] / LimitX;
	screen[1] = obj.element[0][2] / LimitY;
}

}