#include"raytracer.h"
int main() {
	Raytracer* raytracer = new Raytracer;
	raytracer->SetInput( "myscene.txt" );
	raytracer->SetOutput( "mypicture.bmp" );
	//raytracer->Run();
	raytracer->MultiThreadRun();
	//raytracer->DebugRun(740,760,410,430);
	return 0;
}
