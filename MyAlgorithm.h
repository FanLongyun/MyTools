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
	DCM.element[0][1] = sinf(rad_phi) * sinf(rad_theta) * cosf(rad_psi) - sinf(rad_psi) * cosf(rad_phi);
	DCM.element[0][2] = cosf(rad_psi) * sinf(rad_theta) * cosf(rad_phi) + sinf(rad_psi) * sinf(rad_phi);
	DCM.element[1][0] = sinf(rad_phi) * cosf(rad_theta);
	//DCM.element[1][0] = sinf(rad_psi) * cosf(rad_theta);
	DCM.element[1][1] = sinf(rad_psi) * sinf(rad_theta) * sinf(rad_phi) + cosf(rad_psi) * cosf(rad_phi);
	DCM.element[1][2] = sinf(rad_psi) * sinf(rad_theta) * cosf(rad_phi) - cosf(rad_psi) * sinf(rad_phi);
	DCM.element[2][0] = -1.0f * sinf(rad_theta);
	DCM.element[2][1] = cosf(rad_theta) * sinf(rad_phi);
	DCM.element[2][2] = cosf(rad_theta) * cosf(rad_phi);

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
    void normalize();							// 向量标准化
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

}