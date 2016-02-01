#include <sdk.h> // Code::Blocks SDK
#include <configurationpanel.h>
#include "infographie.h"
#include "tgaimage.h"
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <string.h>
// Register the plugin with Code::Blocks.
// We are using an anonymous namespace so we don't litter the global one.
namespace
{
    PluginRegistrant<infographie> reg(_T("infographie"));
}

using namespace std;

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
  bool steep = std::abs(y1-y0) > std::abs(x1-x0) ;
    if(steep){
      std::swap(x0,y0);
      std::swap(x1,y1);
    }
  if(x1<x0){
    std::swap(x0,x1);
    std::swap(y0,y1);
  }
  int dx = x1-x0;
  int dy = y1-y0;
  int derror = std::abs(dy)*2;
  int error = 0;
  int y = y0;
  for (int x = x0; x<x1; x++){
  if(steep){
     image.set(y, x, color);
  }else{
     image.set(x, y, color);
  }
  error += 2*derror;
  if(error > dx){
    y += (y1>y0?1:-1);
    error -= 2*dx;
  }
  }
}


// constructor
infographie::infographie()
{
    // Make sure our resources are available.
    // In the generated boilerplate code we have no resources but when
    // we add some, it will be nice that this code is in place already ;)
    if(!Manager::LoadResource(_T("infographie.zip")))
    {
        NotifyMissingFile(_T("infographie.zip"));
    }
}

// destructor
infographie::~infographie()
{
}

void infographie::OnAttach()
{
    // do whatever initialization you need for your plugin
    // NOTE: after this function, the inherited member variable
    // m_IsAttached will be TRUE...
    // You should check for it in other functions, because if it
    // is FALSE, it means that the application did *not* "load"
    // (see: does not need) this plugin...
}

void infographie::OnRelease(bool appShutDown)
{
    // do de-initialization for your plugin
    // if appShutDown is true, the plugin is unloaded because Code::Blocks is being shut down,
    // which means you must not use any of the SDK Managers
    // NOTE: after this function, the inherited member variable
    // m_IsAttached will be FALSE...
}

int infographie::Execute()
{
    // do your magic ;)
    NotImplemented(_T("infographie::Execute()"));
    return -1;
}

int main(int argc, char** argv) {
	TGAImage image(800, 800, TGAImage::RGB);
	std::ifstream fichier("obj/african_head.obj",std::ios::in);
	if(fichier){
	  printf("Lecture du fichier : \n");
	  string ligne;
	  string x[2942],y[2942],z[2942];
      size_t str,str2,str3;
	  int i=0;
	  int j=0;
	  double x1[1258],y1[1258],z1[1258];

	  while( getline(fichier, ligne)){
	    if(ligne.size()>0){
	      char id= ligne.at(0);
	      char f='f';
	      char v='v';
	      if(id==v){
		cout << ligne << endl;
		/*
		ligne  >> x[j] >>y[j]>>z[j] ;

		cout<<x[j] << endl;
		cout<< y[j]<<endl;
		cout<<z[j]<<endl;
		*/
		/*
		str=ligne.find(' ');
		str2=ligne.find(' ',str+1);
		x1[j]=std::stod(ligne.substr(str,str2));
		str3=ligne.find(' ',str2+1);
		y1[j]=ligne.substr(str2,str3);
		z1[j]=ligne.substr(str3,ligne.length());
		*/
		j++;
	      }
	      if(id==f){
		cout << ligne << endl;
		str=ligne.find(' ');
		//cout << str << endl;

		str2=ligne.find(' ', str+1);
		//cout << str2 << endl ;
		x[i]=ligne.substr(str,str2);
		//cout << x << endl ;
		str3=ligne.find(' ', str2+1);
		y[i]=ligne.substr(str2,str3);
		//cout << y << endl ;
		z[i]=ligne.substr(str3,ligne.length());
		//cout << z << endl;
		str=x[i].find('/');
		x[i]=x[i].substr(0,str);
		str=y[i].find('/');
		y[i]=y[i].substr(0,str);
		str=z[i].find('/');
		z[i]=z[i].substr(0,str);
		cout << x[i] << endl ;
		cout << y[i] << endl ;
		cout << z[i] << endl ;
		i++;
	      }
	    }
	  }


	  fichier.close();
	  printf("Fin de lecture du fichier \n");
	}else{

	}


	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("output.tga");
	return 0;
}
