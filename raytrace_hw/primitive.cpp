#include"primitive.h"
#include<sstream>
#include<cstdio>
#include<string>
#include<cmath>
#include<iostream>
#include<cstdlib>
#include <algorithm>
#include <tuple>
#define ran() ( double( rand() % 32768 ) / 32768 )

const int BEZIER_MAX_DEGREE = 5;
const int Combination[BEZIER_MAX_DEGREE + 1][BEZIER_MAX_DEGREE + 1] =
{	0, 0, 0, 0, 0, 0,
	1, 1, 0, 0, 0, 0,
	1, 2, 1, 0, 0, 0,
	1, 3, 3, 1, 0, 0,
	1, 4, 6, 4, 1, 0,
	1, 5, 10,10,5, 1
};

const int MAX_COLLIDE_TIMES = 10;
const int MAX_COLLIDE_RANDS = 10;


std::pair<double, double> ExpBlur::GetXY()
{
	double x,y;
	x = ran();
	x = pow(2, x)-1;
	y = rand();
	return std::pair<double, double>(x*cos(y),x*sin(y));
}

Material::Material() {
	color = absor = Color();
	refl = refr = 0;
	diff = spec = 0;
	rindex = 0;
	drefl = 0;
	texture = NULL;
	blur = new ExpBlur();
}

void Material::Input( std::string var , std::stringstream& fin ) {
	if ( var == "color=" ) color.Input( fin );
	if ( var == "absor=" ) absor.Input( fin );
	if ( var == "refl=" ) fin >> refl;
	if ( var == "refr=" ) fin >> refr;
	if ( var == "diff=" ) fin >> diff;
	if ( var == "spec=" ) fin >> spec;
	if ( var == "drefl=" ) fin >> drefl;
	if ( var == "rindex=" ) fin >> rindex;
	if ( var == "texture=" ) {
		std::string file; fin >> file;
		texture = new Bmp;
		texture->Input( file );
	}
	if ( var == "blur=" ) {
		std::string blurname; fin >> blurname;
		if(blurname == "exp")
			blur = new ExpBlur();
	}
}

Primitive::Primitive() {
	sample = rand();
	material = new Material;
	next = NULL;
}

Primitive::Primitive( const Primitive& primitive ) {
	*this = primitive;
	material = new Material;
	*material = *primitive.material;
}

Primitive::~Primitive() {
	delete material;
}

void Primitive::Input( std::string var , std::stringstream& fin ) {
	material->Input( var , fin );
}

Sphere::Sphere() : Primitive() {
	//vector3表示一个空间向量
	De = Vector3( 0 , 0 , 1 );
	Dc = Vector3( 0 , 1 , 0 );
}

void Sphere::Input( std::string var , std::stringstream& fin ) {
	//球体的圆心
	if ( var == "O=" ) O.Input( fin );
	//球体的半径
	if ( var == "R=" ) fin >> R;
	//?? De, Dc表示什么？
	if ( var == "De=" ) De.Input( fin );
	if ( var == "Dc=" ) Dc.Input( fin );
	Primitive::Input( var , fin );
}

/*
ray_O: 光线源点
ray_V: 光线的方向
*/
CollidePrimitive Sphere::Collide( Vector3 ray_O , Vector3 ray_V ) {
	//得到单位化的光线方向向量
	ray_V = ray_V.GetUnitVector();
	//注意P的方向，
	Vector3 P = ray_O - O;
	//向量的点乘，得到p在光线方向向量上的投影长度
	double b = -P.Dot( ray_V );
	 
	//
	double det = b * b - P.Module2() + R * R;
	CollidePrimitive ret;

	//det>0时不是一定相交了吗？
	if ( det > EPS ) {
		//得到物体内相交光线的单位长度t
		det = sqrt( det );
		//x1为相交最近距离，x2为相交最远距离
		double x1 = b - det  , x2 = b + det;

		//这一种情况是？？
		if ( x2 < EPS ) return ret;//返回默认值

		//取相交最近点
		if ( x1 > EPS ) {
			//光源点到球体相交点的距离
			ret.dist = x1;
			ret.front = true;
		} 
		//取相交最远点
		else {
			ret.dist = x2;
			ret.front = false;
		} 
	} else 
		return ret;

	//相交点
	ret.C = ray_O + ray_V * ret.dist;
	//相交点的法向量
	ret.N = ( ret.C - O ).GetUnitVector();
	//调整法向量的方向
	if ( ret.front == false ) ret.N = -ret.N;
	ret.isCollide = true;
	ret.collide_primitive = this;
	return ret;
}

//？                                                                                                                                                             
Color Sphere::GetTexture(Vector3 crash_C) {
	Vector3 I = ( crash_C - O ).GetUnitVector();
	//the angle of y axis
	double a = acos( -I.Dot( De ) );
	//the angle of x axis
	double b = acos( std::min( std::max( I.Dot( Dc ) / sin( a ) , -1.0 ) , 1.0 ) );
	double u = a / PI , v = b / 2 / PI;
	//turnover the v coordinate 
	if ( I.Dot( Dc * De ) < 0 ) v = 1 - v;
	return material->texture->GetSmoothColor( u , v );
}


void Plane::Input( std::string var , std::stringstream& fin ) {
	if ( var == "N=" ) N.Input( fin );
	if ( var == "R=" ) fin >> R;
	if ( var == "Dx=" ) Dx.Input( fin );
	if ( var == "Dy=" ) Dy.Input( fin );
	Primitive::Input( var , fin );
}

CollidePrimitive Plane::Collide( Vector3 ray_O , Vector3 ray_V ) {
	ray_V = ray_V.GetUnitVector();
	//将平面法向量单位化？
	N = N.GetUnitVector();
	//得到光线方向向量在法向量方向上的单位投影长度
	double d = N.Dot( ray_V );
	CollidePrimitive ret;
	//
	if ( fabs( d ) < EPS ) return ret;

	//
	double l = ( N * R - ray_O ).Dot( N ) / d;
	//光线与平面没有相交
	if ( l < EPS ) return ret;

	//光线与平面相交
	ret.dist = l;

	//单位投影d<0, 光线射入平面；否则光线射出平面
	ret.front = ( d < 0 );
	//相交光线值
	ret.C = ray_O + ray_V * ret.dist;

	//注意相交光线的法向量方向
	ret.N = ( ret.front ) ? N : -N;
	ret.isCollide = true;
	ret.collide_primitive = this;
	return ret;
}

Color Plane::GetTexture(Vector3 crash_C) {
	double u = crash_C.Dot( Dx ) / Dx.Module2();
	double v = crash_C.Dot( Dy ) / Dy.Module2();
	return material->texture->GetSmoothColor( u , v );
}

void Square::Input( std::string var , std::stringstream& fin ) {
	if ( var == "O=" ) O.Input( fin );
	//Dx为正方体x轴的向量
	if ( var == "Dx=" ) Dx.Input( fin );
	//Dy为正方体y轴的向量
	if ( var == "Dy=" ) Dy.Input( fin );
	Primitive::Input( var , fin );
}

//(ray_o+ray_v*t - O).N=0
//N = Dx x Dy
CollidePrimitive Square::Collide( Vector3 ray_O , Vector3 ray_V ) {
	CollidePrimitive ret;
	ray_V = ray_V.GetUnitVector();
	/*Vector3 N = Dx * Dy;
	Vector3 P = ray_O - O;
	double t = (-P.Dot(N) / (ray_V.Dot(N)));

	if (fabs(t) < EPS)
		return ret;
	Vector3 iscollideP = ray_O + ray_V * t;
	Vector3 V = iscollideP - O;

	double projx = V.Dot(Dx);
	double projy = V.Dot(Dy);
	if (projx < Dx.Module2() && projx > 0 &&
		projy < Dy.Module2() && projy > 0)
	{
		ret.C = iscollideP;
		ret.front = (ray_V.Dot(N) < 0);
		ret.dist = t;
		ret.N = (ret.front) ? N : -N;
		ret.isCollide = true;
		ret.collide_primitive = this;
	}
	return ret;*/

	//正方体边长
	double len = Dx.Module();
	//正方体上最小点
	Vector3 minP = Vector3(O.x - len / 2, O.y - len / 2, O.z - len / 2);
	//正方体上最大点
	Vector3 maxP = Vector3(O.x + len / 2, O.y + len / 2, O.z + len / 2);
	Vector3 isCollideP;

	//Dx cross product Dy
	Vector3 Dz = Dx * Dy;
	Vector3 N;
	
	ray_V = ray_V.GetUnitVector();
	double t = 0;

	//分别判断光线与各面的相交情况
	//使用射线与平面相交公式计算交点
	/*
		ray = ray_O + ray_V * t
		point * N = D [point is on the plane, N is the plane's normal, D is the dx,dy or dz]
		t = (D - ray_O * N) / (ray_V * N)
		t = (Dx -ray_O.x) / ray_V.x; (y,z as the same)
	*/

	//光线是否与x轴平行
	if (ray_V.x != 0)
	{

		if (ray_V.x > 0)//光线沿x轴正方向偏移
		{
			t = (minP.x - ray_O.x) / ray_V.x;
			N = -Dx.GetUnitVector();
		}


		else//光线沿x轴负方向偏移
		{
			t = (maxP.x - ray_O.x) / ray_V.x;
			N = Dx.GetUnitVector();
		}


		if (t > 0)
		{
			//交点坐标
			isCollideP = ray_O + ray_V * t;

			//判断交点是否在正方体内
			if (minP.y < isCollideP.y && isCollideP.y < maxP.y &&
				minP.z < isCollideP.z && isCollideP.z < maxP.z)
			{
				ret.dist = t;
				ret.C = isCollideP;

				//front如何判断？
				ret.front = (ray_V.Dot(N) < 0);
				ret.N = (ret.front) ? N : -N;
				ret.isCollide = true;
				ret.collide_primitive = this;
				return ret;
			}
	
		}

	}

	if (ray_V.y != 0.0)
	{
		if (ray_V.y > 0)
		{
			t = (minP.y - ray_O.y) / ray_V.y;
			N = -Dy.GetUnitVector();
		}

		else
		{
			t = (maxP.y - ray_O.y) / ray_V.y;
			N = Dy.GetUnitVector();
		}

		if (t > 0.0)
		{
			isCollideP = ray_O + ray_V * t;

			if (minP.z < isCollideP.z && isCollideP.z < maxP.z &&
				minP.x < isCollideP.x && isCollideP.x < maxP.x)
			{
				ret.dist = t;
				ret.C = isCollideP;
				
				//front?
				ret.front = (ray_V.Dot(N) < 0);
				ret.N = (ret.front) ? N : -N;
				ret.isCollide = true;
				ret.collide_primitive = this;
				return ret;
			}

		}
	}

	if (ray_V.z != 0.0)
	{
		if (ray_V.z > 0.0)
		{
			t = (minP.z - ray_O.z) / ray_V.z;
			N = -Dz.GetUnitVector();
		}

		else
		{
			t = (maxP.z - ray_O.z) / ray_V.z;
			N = Dz.GetUnitVector();
		}

		if (t > 0.0)
		{
			isCollideP = ray_O + ray_V * t;

			if (minP.x < isCollideP.x && isCollideP.x < maxP.x &&
				minP.y < isCollideP.y && isCollideP.y < maxP.y)
			{
				ret.dist = t;
				ret.C = isCollideP;

				//front?
				ret.front = (ray_V.Dot(N) < 0);
				ret.N = (ret.front) ? N : -N;
				ret.isCollide = true;
				ret.collide_primitive = this;
				return ret;
			}
		}

	}
	return ret;
}

Color Square::GetTexture(Vector3 crash_C) {
	double u = (crash_C - O).Dot( Dx ) / Dx.Module2() / 2 + 0.5;
	double v = (crash_C - O).Dot( Dy ) / Dy.Module2() / 2 + 0.5;
	return material->texture->GetSmoothColor( u , v );
}

void Cylinder::Input( std::string var , std::stringstream& fin ) {
	if ( var == "O1=" ) O1.Input( fin );
	if ( var == "O2=" ) O2.Input( fin );
	if ( var == "R=" ) fin>>R; 
	Primitive::Input( var , fin );
}

//calculate the cylinder's top and bottom intersect point
static void intersect(Vector3 ray_O, Vector3 ray_V, Vector3 N,
	Vector3 p,double &distance, bool &front, Vector3 &intersect)
{
	double dist = (ray_O - p).Dot(N.GetUnitVector());
	front = dist >= 0;
	distance = dist / (-ray_V.GetUnitVector().Dot(N.GetUnitVector()));
	intersect = ray_O + ray_V.GetUnitVector() * distance;
}

CollidePrimitive Cylinder::Collide( Vector3 ray_O , Vector3 ray_V ) {

	CollidePrimitive ret;
	ray_V = ray_V.GetUnitVector();

	//use tuple <distance, front, intersect>
	using collide_tuple = std::tuple<double, bool, Vector3>;
	collide_tuple o1, o2, c1, c2;

	//cylinder's bottom and top

	//calculate the O1's plane intersection
	Vector3 axis = O2 - O1;//the cylinder's height axis
	intersect(ray_O, ray_V, axis, O1, std::get<0>(o1),
		std::get<1>(o1), std::get<2>(o1));

	//whether the intersection is within the O1's plane
	if ((std::get<2>(o1) - O1).Module2() >= R * R)
		//distance
		std::get<0>(o1) = -1;
	
	//O2's plane
	intersect(ray_O, ray_V, axis, O2, std::get<0>(o2),
		std::get<1>(o2), std::get<2>(o2));
	if ((std::get<2>(o2) - O2).Module2() >= R * R)
		std::get<0>(o2) = -1;

	//calculate the cylinder's side plane intersection
	/*
	cycle:|point-O|=R
	|(o+vt-O1) x N| = R
	|(o-O1) x N + vtN| = R
	[x,y,z]=(o-O1) x N; [p,q,r]=vtN
	*/

	Vector3 N = axis.GetUnitVector();
	Vector3 pn = (ray_O - O1) * N;
	Vector3 vn = ray_V * N;

	double x = pn.x, y = pn.y, z = pn.z;
	double p = vn.x, q = vn.y, r = vn.z;
	double a = p * p + q * q + r * r;
	double b = 2 * (p * x + q * y + r * z);
	double c = x * x + y * y + z * z - R * R;

	double delta = b * b - 4 * a * c;
	if (delta < EPS)
		return ret;
	delta = sqrt(delta);
	double t1 = (-b - delta) / 2 / a;
	double t2 = (-b + delta) / 2 / a;
	if (t1 > EPS)
	{
		std::get<0>(c1) = t1 * ray_V.Module();
		std::get<2>(c1) = ray_O + ray_V * t1;
	}
	if (t2 > EPS)
	{
		std::get<0>(c2) = t2 * ray_V.Module();
		std::get<2>(c2) = ray_O + ray_V * t2;
	}
	//the camera is outside the cylinder
	if (t1 > EPS && t2 > EPS)
	{
		std::get<1>(c1) = (t1 <= t2);
		std::get<1>(c2) = (t2 <= t1);
	}
	//camera is outside the cylinder
	else
	{
		std::get<1>(c1) = std::get<1>(c2) = false;
	}

	//whether the point is within the cylinder
	if (std::get<0>(c1) >= EPS)
	{
		double h = (std::get<2>(c1) - O1).Dot(axis);
		if (h<0 || h>axis.Module2())
			std::get<0>(c1) = -1;
	}
	if (std::get<0>(c2) >= EPS)
	{
		double h = (std::get<2>(c2) - O1).Dot(axis);
		if (h<0 || h>axis.Module2())
			std::get<0>(c2) = -1;
	}

	//choose the smallest distance (dist>0)
	collide_tuple* mini = std::min({ &o1,&o2,&c1,&c2 },
		[](const collide_tuple* a, const collide_tuple* b)
		{
			//return (a<b)?a:b
			if (std::get<0>(*a) < EPS)
				return false;
			if (std::get<0>(*b) < EPS)
				return true;
			return std::get<0>(*a) < std::get<0>(*b);
		});

	//the point is outside the cylinder
	if (std::get<0>(*mini) < EPS)
		return ret;

	ret.isCollide = true;
	ret.dist = std::get<0>(*mini);
	ret.front = std::get<1>(*mini);
	ret.C = std::get<2>(*mini);
	if (mini == &o1)
		ret.N = (O1 - O2).GetUnitVector();
	else if (mini == &o2)
		ret.N = (O2 - O1).GetUnitVector();
	else
		ret.N = (axis * ((std::get<2>(*mini) - O1) * axis)).GetUnitVector();
	if (ret.front == false)
		ret.N = -ret.N;

	ret.collide_primitive = this;
	return ret;

}

Color Cylinder::GetTexture(Vector3 crash_C) {
	Vector3 O = (O1 + O2) / 2;
	Vector3 De = (O2 - O1).GetUnitVector();
	Vector3 Dc = (O2 - O1).GetAnVerticalVector().GetUnitVector();
	Vector3 I = (crash_C - O).GetUnitVector();
	double a = acos(-I.Dot(De));
	double b = acos(std::min(std::max(I.Dot(Dc) / sin(a), -1.0), 1.0));
	double u = a / PI, v = b / 2 / PI;
	if (I.Dot(Dc * De) < 0)
		v = 1 - v;
	return material->texture->GetSmoothColor( u , v );
}

void Bezier::Input( std::string var , std::stringstream& fin ) {
	if ( var == "O1=" ) O1.Input( fin );
	if ( var == "O2=" ) O2.Input( fin );
	if ( var == "P=" ) {
		degree++;
		double newR, newZ;
		fin>>newZ>>newR;
		R.push_back(newR);
		Z.push_back(newZ);
	}
	if ( var == "Cylinder" ) {
		double maxR = 0;
		for(int i=0;i<R.size();i++)
			if(R[i] > maxR)
				maxR = R[i];
		boundingCylinder = new Cylinder(O1, O2, maxR);
		N = (O1 - O2).GetUnitVector();
		Nx = N.GetAnVerticalVector();
		Ny = N * Nx;
	}
	Primitive::Input( var , fin );
}

CollidePrimitive Bezier::Collide( Vector3 ray_O , Vector3 ray_V ) {
	CollidePrimitive ret;
	//NEED TO IMPLEMENT
	return ret;
}

Color Bezier::GetTexture(Vector3 crash_C) {
	double u = 0.5 ,v = 0.5;
	//NEED TO IMPLEMENT
	return material->texture->GetSmoothColor( u , v );
}

