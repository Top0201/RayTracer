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
	double refl , refr;//����ϵ��������ϵ��
	double diff , spec;//������ϵ�������淴��ϵ��
	double rindex;
	double drefl;//��������ʲô��
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
	//�õ�����׷�ٺ����ɫ
	virtual Color GetTexture(Vector3 crash_C) = 0;
	virtual bool IsLightPrimitive(){return false;}
};

struct CollidePrimitive {
	//�Ƿ�����ײ
	bool isCollide;
	//������ײ����������
	Primitive* collide_primitive;
	//N��������C�����ڵ���ײ��������
	Vector3 N , C;

	//distance��
	double dist;
	//�Ƿ�ǰ�棿
	bool front;
	CollidePrimitive(){isCollide = false; collide_primitive = NULL; dist = BIG_DIST;}

	//�õ���ײ�����ɫ
	Color GetTexture(){return collide_primitive->GetTexture(C);}
};

//sphere is a child for primitive
class Sphere : public Primitive {
protected:
	//De, Dc����������ɫ���������ʲô��
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
	//N��ƽ��ķ�����
	//Dx��Dy�����������׷�ٵ���ɫ
	Vector3 N , Dx , Dy;
	//ƽ��İ뾶��Χ��
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
	//OΪ����ԭ�㣬Dx��Dy��
	//Dx,DyΪ�������x���y��
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
	//Բ���������Բ��
	Vector3 O1, O2;
	//�뾶
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
