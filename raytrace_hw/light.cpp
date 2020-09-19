#include"light.h"
#include<sstream>
#include<string>
#include<cmath>
#include<cstdlib>
#define ran() ( double( rand() % 32768 ) / 32768 )

Light::Light() {
	sample = rand();
	next = NULL;
	lightPrimitive = NULL;
}

void Light::Input( std::string var , std::stringstream& fin ) {
	if ( var == "color=" ) color.Input( fin );
}

void PointLight::Input( std::string var , std::stringstream& fin ) {
	if ( var == "O=" ) O.Input( fin );
	Light::Input( var , fin );
}


double PointLight::CalnShade( Vector3 C , Primitive* primitive_head , int shade_quality ) {
	Vector3 V = O - C;
	double dist = V.Module();
	//���������ޱ��ڵ�������shadowϵ��
	for ( Primitive* now = primitive_head ; now != NULL ; now = now->GetNext() )
	{
		CollidePrimitive tmp = now->Collide(C, V);
		if ( EPS < (dist - tmp.dist) ) 
			//�ڵ�shadow=0
			return 0;
	}
	//δ�ڵ�shadow=1
	return 1;
}

void SquareLight::Input( std::string var , std::stringstream& fin ) {
	if ( var == "O=" ) O.Input( fin );
	if ( var == "Dx=" ) Dx.Input( fin );
	if ( var == "Dy=" ) Dy.Input( fin );
	Light::Input( var , fin );
}

//?���ʵ��area light
double SquareLight::CalnShade( Vector3 C , Primitive* primitive_head , int shade_quality ) {
	return 1;
	
}

Primitive* SquareLight::CreateLightPrimitive()
{
	//��OΪ������ƽ�����ĵ�ƽ���
	PlaneAreaLightPrimitive* res = new PlaneAreaLightPrimitive(O, Dx, Dy, color);
	lightPrimitive = res;
	return res;
}



void SphereLight::Input( std::string var , std::stringstream& fin ) {
	if ( var == "O=" ) O.Input( fin );
	if ( var == "R=" ) fin>>R;
	Light::Input( var , fin );
}

//�㷨��������θĽ���
double SphereLight::CalnShade( Vector3 C , Primitive* primitive_head , int shade_quality ) {
	long long shade = 0;
	//NEED TO IMPLEMENT
	double stackAngle = PI / 3;
	double sectorAngle = 2 * PI / 3;

	long long lightCount = 0;
	long long shadowCount = 0;
	//�����������
	for (int i = 0; i <= PI; i += stackAngle)
	{
		for (int j = 0; j <= 2 * PI; j += sectorAngle)
		{
			double x = R * sin(i) * sin(j);
			double y = R * sin(i) * cos(j);
			double z = R * cos(i);
			Vector3 P = Vector3(x, y, z);
			Vector3 V = P - C;
			double dist = V.Module();

			for (Primitive* now = primitive_head; now != NULL; now = now->GetNext())
			{
				CollidePrimitive tmp = now->Collide(C, V);
				if (EPS < (dist - tmp.dist))
					shadowCount += 1;
				else
					lightCount += 1;
			}

		}
	}

	shade = lightCount / (lightCount + shadowCount);
	return shade;
}


Primitive* SphereLight::CreateLightPrimitive()
{
	SphereLightPrimitive* res = new SphereLightPrimitive(O, R, color);
	lightPrimitive = res;
	return res;
}

