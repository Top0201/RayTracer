#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include"color.h"
#include"vector3.h"
#include"bmp.h"
#include<iostream>
#include<sstream>
#include<string>
#include<vector>

extern const double EPS;
extern const double PI;
const double BIG_DIST = 1e100;

class Blur {
public:
	virtual std::pair<double, double> GetXY() = 0;
};

class ExpBlur : public Blur {
public:
	std::pair<double, double> GetXY();
};

class Material {
public:
	Color color , absor;
	double refl , refr;//反射系数，折射系数
	double diff , spec;//漫反射系数，镜面反射系数
	double rindex;
	double drefl;//？？这是什么？
	Bmp* texture;
	Blur* blur;

	Material();
	~Material() {}

	void Input( std::string , std::stringstream& );
};

struct CollidePrimitive;

class Primitive {
protected:
	int sample;
	Material* material;
	Primitive* next;

public:

	Primitive();
	Primitive( const Primitive& );
	~Primitive();
	
	int GetSample() { return sample; }
	Material* GetMaterial() { return material; }
	Primitive* GetNext() { return next; }
	void SetNext( Primitive* primitive ) { next = primitive; }

	virtual void Input( std::string , std::stringstream& );
	virtual CollidePrimitive Collide( Vector3 ray_O , Vector3 ray_V ) = 0;
	//得到光线追踪后的颜色
	virtual Color GetTexture(Vector3 crash_C) = 0;
	virtual bool IsLightPrimitive(){return false;}
};

struct CollidePrimitive {
	//是否发生碰撞
	bool isCollide;
	//发生碰撞的物体类型
	Primitive* collide_primitive;
	//N法向量，C物体内的碰撞光线向量
	Vector3 N , C;

	//distance？
	double dist;
	//是否前面？
	bool front;
	CollidePrimitive(){isCollide = false; collide_primitive = NULL; dist = BIG_DIST;}

	//得到碰撞后的颜色
	Color GetTexture(){return collide_primitive->GetTexture(C);}
};

//sphere is a child for primitive
class Sphere : public Primitive {
protected:
	//De, Dc用来计算颜色的输出，是什么？
	Vector3 O , De , Dc;
	double R;

public:
	Sphere();
	~Sphere() {}

	void Input( std::string , std::stringstream& );
	CollidePrimitive Collide( Vector3 ray_O , Vector3 ray_V );
	Color GetTexture(Vector3 crash_C);
};

class SphereLightPrimitive : public Sphere{
public:
	SphereLightPrimitive(Vector3 pO, double pR, Color color) : Sphere()
	{O = pO; R = pR; material->color = color; }
	bool IsLightPrimitive(){return true;}
};

class Plane : public Primitive {
protected:
	//N：平面的法向量
	//Dx，Dy用来计算光线追踪的颜色
	Vector3 N , Dx , Dy;
	//平面的半径范围？
	double R;

public:
	Plane() : Primitive() {}
	~Plane() {}

	void Input( std::string , std::stringstream& );
	CollidePrimitive Collide( Vector3 ray_O , Vector3 ray_V );
	Color GetTexture(Vector3 crash_C);
};

class Square : public Primitive {
protected:
	//O为坐标原点，Dx，Dy？
	//Dx,Dy为正方体的x轴和y轴
	Vector3 O , Dx , Dy;

public:
	Square() : Primitive() {}
	~Square() {}

	void Input( std::string , std::stringstream& );
	CollidePrimitive Collide( Vector3 ray_O , Vector3 ray_V );
	Color GetTexture(Vector3 crash_C);
};

class PlaneAreaLightPrimitive : public Square{
public:
	PlaneAreaLightPrimitive(Vector3 pO, Vector3 pDx, Vector3 pDy, Color color): Square()
	{O = pO; Dx = pDx; Dy = pDy; material->color = color; }
	bool IsLightPrimitive(){return true;}
};

class Cylinder : public Primitive {
	//圆柱上下面的圆心
	Vector3 O1, O2;
	//半径
	double R;

public:
	Cylinder() : Primitive() {}
	Cylinder(Vector3 pO1, Vector3 pO2, double pR) : Primitive() {O1 = pO1; O2 = pO2; R = pR; }
	~Cylinder() {}

	void Input( std::string , std::stringstream& );
	CollidePrimitive Collide( Vector3 ray_O , Vector3 ray_V );
	Color GetTexture(Vector3 crash_C);
};

class Bezier : public Primitive {
	Vector3 O1, O2;
	Vector3 N, Nx, Ny;
	std::vector<double> R;
	std::vector<double> Z;
	int degree;
	Cylinder* boundingCylinder;

public:
	Bezier() : Primitive() {boundingCylinder = NULL; degree = -1;}
	~Bezier() {}

	void Input( std::string , std::stringstream& );
	CollidePrimitive Collide( Vector3 ray_O , Vector3 ray_V );
	Color GetTexture(Vector3 crash_C);
};

#endif
