
/*-------------------------------------------------------------------------*\

This file is part of libvoxel, a C++ template library for handelling 3D images.

Developed by:
 - Ali Q Raeini (2010-2022)

You can redistribute this code and/or modify this code under the
terms of the GNU General Public License (GPL) as published by the
Free Software Foundation, either version 3 of the License, or (at
your option) any later version. see <http://www.gnu.org/licenses/>.

\*-------------------------------------------------------------------------*/

#include <math.h>
#include <fstream>
#include <assert.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include "voxelImage.h"

#include <functional>   // std::mem_fn

								namespace MCTProcessing _begins_
using namespace std;

struct shape {
 public:
	constexpr static int invalidv = 0x80000000; // 2147483648 largest negative
	const char polyType; // shape // TODO use enum instead (?)
	int insidev=0;

	shape(char _polyType): polyType(_polyType), insidev(0){}
	shape(char _polyType, int _insidev): polyType(_polyType), insidev(_insidev){}

	int isBefore(dbl3) const { return false; };
	int isAfter(dbl3) const { return false; };

	template<typename T, class Shape>
	static void setIn(voxelImageT<T>& vImg, const Shape& shp) {
		(cout <<__FUNCTION__<<'_'<<shp.polyType<<": "<<shp.insidev<<",  ").flush();
		dbl3 xmin = vImg.X0(),  dx = vImg.dx();
		forAllkji_(vImg) {
			int vv = shp.value(dbl3(i+0.5,j+0.5,k+0.5)*dx+xmin);
			if(vv!=shape::invalidv) vImg(i,j,k)  = vv;
		}
	}
	template<typename T, class Shape>
	static void addTo(voxelImageT<T>& vImg, const Shape& shp) {
		(cout <<__FUNCTION__<<'_'<<shp.polyType<<": "<<shp.insidev<<",  ").flush();
		dbl3 xmin = vImg.X0(),  dx = vImg.dx();
		forAllkji_(vImg)
			if(int vv = shp.value(dbl3(i+0.5,j+0.5,k+0.5)*dx+xmin); vv!=shape::invalidv) vImg(i,j,k) += vv;
	}
	template<typename T, class Shape>
	static void setBefor(voxelImageT<T>& vImg, const Shape& shp) {
		(cout <<__FUNCTION__<<'_'<<shp.polyType<<": "<<shp.insidev<<",  ").flush();
		dbl3 xmin = vImg.X0(),  dx = vImg.dx();
		forAllkji_(vImg)
			if(shp.isBefore(dbl3(i+0.5,j+0.5,k+0.5)*dx+xmin))  vImg(i,j,k)  = shp.insidev;
	}
	template<typename T, class Shape>
	static void setAfter(voxelImageT<T>& vImg, const Shape& shp) {
		(cout <<__FUNCTION__<<'_'<<shp.polyType<<": "<<shp.insidev<<",  ").flush();
		dbl3 xmin = vImg.X0(),  dx = vImg.dx();
		forAllkji_(vImg)
			if(shp.isAfter (dbl3(i+0.5,j+0.5,k+0.5)*dx+xmin))  vImg(i,j,k) = shp.insidev;
	}
	template<typename T, class Shape>
	static void addBefor(voxelImageT<T>& vImg, const Shape& shp) {
		(cout <<__FUNCTION__<<'_'<<shp.polyType<<": "<<shp.insidev<<",  ").flush();
		dbl3 xmin = vImg.X0(),  dx = vImg.dx();
		forAllkji_(vImg)
			if(shp.isBefore(dbl3(i+0.5,j+0.5,k+0.5)*dx+xmin))  vImg(i,j,k) += shp.insidev;
	}
	template<typename T, class Shape>
	static void addAfter(voxelImageT<T>& vImg, const Shape& shp) {
		(cout <<__FUNCTION__<<'_'<<shp.polyType<<": "<<shp.insidev<<",  ").flush();
		dbl3 xmin = vImg.X0(),  dx = vImg.dx();
		forAllkji_(vImg)
			if(shp.isAfter(dbl3(i+0.5,j+0.5,k+0.5)*dx+xmin))  vImg(i,j,k) += shp.insidev;
	}
};


class cylinder : public shape {
	dbl3 p1, p2;  double rr;
	double mag_p12;
 public:
	cylinder(stringstream & ins) : shape(/*polyType*/'c') {
		//http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
		ins>>p1 >> p2 >>rr >>insidev;
		mag_p12=mag(p2-p1);
		cout <<"cylinder  p1="<<p1<<",   p2="<<p2<<",   r="<<rr<<"   value="<<insidev;
	};

	int value(dbl3 ij) const { return ( mag((ij-p1)^(ij-p2)) <= rr*mag_p12 )  ?  insidev : shape::invalidv; }
 };



class layer : public shape {
	double p1, p2;  //!< p1: location on bottom plane. p2: location on top plane.
	dbl3 nrm; //!< normal vector to planes, along which locations are recorded
  public:
	layer(stringstream & ins): shape(/*polyType*/'l'), nrm(0,1,0) {
		dbl3 Pl; //point through bottom plate,
		ins>>Pl>>nrm>>insidev; p2=mag(nrm); nrm/=p2;  p1=Pl&nrm; p2+=p1;
		cout <<"layer  p1="<<Pl<<",  height="<<p2-p1<<",  normal="<<nrm<<",  insidev="<<insidev;
	}
	int value(dbl3 ij)    const {  return  (ij&nrm)>=p1 && (ij&nrm)<p2  ? insidev : shape::invalidv;  }
	int isBefore(dbl3 ij) const {  return  (ij&nrm) < p1; }// used for network cut
	int isAfter (dbl3 ij) const {  return  (ij&nrm) >= p2;  }
};



 // outdated
class paraPlates : public shape {
	/// input mX dy mZ y1, mX and mZ are slopes, dy is separation, y1 is elevation of bottom plane
	dbl3 p1, p2;

  public:
	paraPlates(stringstream & ins): shape(/*polyType*/'f') {
		ins>>p1;   p2 = p1;   p1.y=0.; ins>>p1.y;
		cout <<"paraPlates slope1,mX="<<p1.x<<"    separation,dy="<<p2.y<<"   slope2,mZ="<<p1.z<<"   shift, y1:"<<p1.y;
	};
	int value(dbl3 ij)   const  {  return ( ij.y > p1.x*ij.x+p1.z*ij.z+p1.y
	                                     && ij.y < p2.x*ij.x+p2.z*ij.z+p2.y )  ? insidev : shape::invalidv;  }
	int isBefore(dbl3 ij) const {  return ( ij.y < p2.x*ij.x+p2.z*ij.z+p2.y ); }// used for network cut
	int isAfter (dbl3 ij) const {  return ( ij.y > p1.x*ij.x+p1.z*ij.z+p1.y );  }
};

class kube : public shape {
	/// TODO  generalize to hexahedron, TODO optimize
	dbl3 p1,p2;
 public:
	kube(stringstream & ins): shape(/*polyType*/'k') {
		ins>>p1>>p2>>insidev;
		cout <<"kube "<<p1<<" + "<<p2<<", value="<<insidev;
		p2+=p1;
	};
	int value(dbl3 ij)  const {  return ij.x >= p1.x && ij.y >= p1.y && ij.z >= p1.z &&   ij.x < p2.x && ij.y < p2.y && ij.z < p2.z
		                                 ? insidev : shape::invalidv;  }
};

class plate : public shape {
	//! plane capped by sphere
	dbl3 p1,p2, po_;
 public:

	plate(stringstream & ins): shape(/*polyType*/'p') {
		cout<<"Error: fix me saqdakjoigfgfgfg "<<endl;
		ins>>p1>>po_; p2=p1;
		cout <<"\n plate: slope_x="<<p1.x<<"  r_cap="<<p1.y<<"  slope_z="<<p1.z<<"\n"<<endl;
		p1.y=0.;    po_.z = po_.y;  po_.y = (p1.x*po_.x+p1.z*po_.z+p1.y);
		cout <<"plate a1="<<p1.x<<"   a2="<<p2.x<<"    r="<<p1.y<<"    b2="<<p2.y<<';';
	};
	int value(dbl3 ij) const {  return  ( ij.y > p1.x*ij.x + p1.z*ij.z+p1.y && mag(ij-po_) <= p2.y )  ?  insidev : shape::invalidv;  }
};

class sphere : public shape {
	dbl3 p1; double r2;
  public:
	sphere(stringstream & ins): shape(/*polyType*/'s') {
		double rr;
		ins>>p1 >>rr >>insidev;   r2=rr*rr;
		cout <<"sphere  c="<<p1<<",  r="<<rr;
	}
	sphere(dbl3 _p1, double rr, int _insidev): shape(/*polyType*/'s',_insidev),
			p1(_p1), r2(rr*rr)
		{}

	int value(dbl3 ij) const {  return  ( magSqr(ij-p1)<r2 )  ?  insidev : shape::invalidv;  }
};



//template<typename T, typename _operate_>
//void applyShapeOper(voxelImageT<T> & vImg, const shape& sh, _operate_& func) {
#define _SHAPERATEPy(vImg, sh, _operate_)  \
	switch (sh.polyType) {                 \
		case 's':  shape::_operate_(vImg, static_cast<const      sphere&>(sh)); break;  \
		case 'p':  shape::_operate_(vImg, static_cast<const       plate&>(sh)); break;  \
		case 'f':  shape::_operate_(vImg, static_cast<const  paraPlates&>(sh)); break;  \
		case 'k':  shape::_operate_(vImg, static_cast<const        kube&>(sh)); break;  \
		case 'l':  shape::_operate_(vImg, static_cast<const       layer&>(sh)); break;  \
		case 'c':  shape::_operate_(vImg, static_cast<const    cylinder&>(sh)); break;  \
		default:  cout <<"\n unregistered shape type: "<<sh.polyType<<"\n"<<endl;        \
	}
//}

#define _SHAPERATE(_operate_)  \
	if(ins.peek()=='?') { ins.str("\"s(phere), p(late, capped), f(flat-plates), c(ylinder) or k(ube)\" position..."); return true; } \
	std::string tmpc;  ins>>tmpc; \
	switch (tmpc[0]) {                 \
		case 's':  shape::_operate_(vImg, sphere(     ins)); break; \
		case 'p':  shape::_operate_(vImg, plate(      ins)); break; \
		case 'f':  shape::_operate_(vImg, paraPlates( ins)); break; \
		case 'k':  shape::_operate_(vImg, kube(       ins)); break; \
		case 'l':  shape::_operate_(vImg, layer(      ins)); break; \
		case 'c':  shape::_operate_(vImg, cylinder(   ins)); break; \
		default:  cout <<"\n unsupported shape type: "<<tmpc<<"\n"<<endl; \
	}

template<typename T>  bool Paint( stringstream& ins, voxelImageT<T> & vImg)  {
	_SHAPERATE(setIn);
	return true;
}

template<typename T>  bool PaintAdd( stringstream& ins, voxelImageT<T> & vImg) {
	_SHAPERATE(addTo);
	return true;
}

template<typename T>  bool PaintBefore( stringstream& ins, voxelImageT<T> & vImg) {
	_SHAPERATE(setBefor);
	return true;
}

template<typename T>  bool PaintAfter( stringstream& ins, voxelImageT<T> & vImg) {
	_SHAPERATE(setAfter);
	return true;
}

template<typename T>  bool PaintAddBefore( stringstream& ins, voxelImageT<T> & vImg) {
	_SHAPERATE(addBefor);
	return true;
}

template<typename T>  bool PaintAddAfter( stringstream& ins, voxelImageT<T> & vImg) {
	_SHAPERATE(addAfter);
	return true;
}

								_end_of_(namespace MCTProcessing)
