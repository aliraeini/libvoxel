#pragma once
/*-------------------------------------------------------------------------*\

This file is part of libvoxel, a C++ template library for handelling 3D images.

Developed by:
 - Ali Q Raeini (2010-2022)

You can redistribute this code and/or modify this code under the
terms of the GNU General Public License (GPL) as published by the
Free Software Foundation, either version 3 of the License, or (at
your option) any later version. see <http://www.gnu.org/licenses/>.

\*-------------------------------------------------------------------------*/


#ifdef _ExtraVxlTypes
#  define SupportedVoxTyps  unsigned char,unsigned short,int,float,short,unsigned int,double,float3,dbl3
#else
#  ifdef _VoxBasic8
#    define SupportedVoxTyps  unsigned char
#  else
#    define SupportedVoxTyps  unsigned char,unsigned short,int,float
#  endif
#endif //_ExtraVxlTypes

#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>
#include <valarray>
#include <cassert>
#include <sstream>
#include <memory>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <functional>

#include "typses.h"

#ifdef _STOR_PUB
#include "SiRun.h" //_STOR
#endif //_STOR_PUB

//! suffix,  set/get default suffix, uses static storage.
inline const std::string& imgExt(const std::string& defSuffix="") {
	#ifdef ZLIB
  	   static std::string defSuffix_=".raw.gz";
	#else
	  #ifdef TIFLIB
	   static std::string defSuffix_=".tif";
	  #else
  	   static std::string defSuffix_=".raw";
	  #endif //TIFLIB
	#endif //ZLIB

	///. set via OutputFormat keyword
	if (defSuffix.size()) {
		if(defSuffix[0]!='.') defSuffix_="."+defSuffix;
		else                  defSuffix_=defSuffix;
		if( defSuffix_!=".tif" && defSuffix_!=".raw.gz" && defSuffix_!=".am" &&
		    defSuffix_!=".raw" && defSuffix_!=".dat" && defSuffix_!=".txt")
			std::cout<<"\nError: wrong default image format: "<<defSuffix_<<"\n"<<std::endl;
	}

	return defSuffix_;
}


template<typename T>
using intOr = typename std::conditional<(sizeof(T) > sizeof(short)),T,int>::type;
#define Tint  intOr<T> // cast to int or larger type

extern int maxNz;

//! 3D voxel data, on Cartesian uniform grids, with efficient access functions, and I/O into  .tif, .raw.gz .am, .raw, .dat file formats.
template <typename T>
class voxelField {
 public:

	long long nij_;
	int3  nnn_;
	std::vector<T> data_;


	voxelField(): nij_(0), nnn_(0,0,0) {}
	voxelField(int3 n)          {  reset(n);  }
	voxelField(int3 n, T value) {  reset(n, value);  }
	voxelField(int n1, int n2, int n3, T value) {  reset(int3(n1,n2,n3), value);  }
	virtual ~voxelField() {}

	void reset(int3 n);
	void reset(int3 n, T value);
	void reset(int n1, int n2, int n3, T value) {  reset(int3(n1,n2,n3), value);  }
	void readMicroCT(std::string);
	bool readAscii(std::string);
	void readAscii(std::ifstream& in);
	int  readBin(std::string fileName, int nSkipBytes=0);
	int  readBin(std::string fileName,int iBgn,int iEndp1 , int jBgn,int jEndp1 , int kBgn,int kEndp1, int nSkipBytes=0);
	void writeNoHdr(std::string fileName) const;
	void writeBin(std::string fileName) const;
	void writeBin(std::string fileName,int iBgn,int iEndp1 , int jBgn,int jEndp1 , int kBgn,int kEndp1 ) const;
	void writeAscii(std::string fileName) const;
	void writeAscii(std::string fileName,int iBgn,int iEndp1 , int jBgn,int jEndp1 , int kBgn,int kEndp1) const;
	void writeRotatedXZ(std::ofstream& of) const;

	void writeHeader(std::string fileName, int3 iBgn, int3 iEnd, dbl3 dx, dbl3 X0) const;

	void setLayer(int k, const T* Values);
	void setSlice(char dir, int ijk, T vv);
	void replacezLayer(int j, int fromj);
	void replaceyLayer(int j, int fromj);
	void replacexLayer(int i, int fromi);
	void setBlock(int n1, int n2, int n3, const voxelField<T>&Values);
	void setFrom(const voxelField<T>&Values, int n1, int n2, int n3);
	void swapData(voxelField<T>& other) noexcept;


	const T& operator()(int i, int j, size_t k) const { return data_[k*nij_+j*nnn_.x+i]; }
	      T& operator()(int i, int j, size_t k)       { return data_[k*nij_+j*nnn_.x+i]; }
	      T& operator()(long long ii, size_t k){ return data_[k*nij_+ii]; }
	const T& operator()(size_t iii)      const { return data_[iii]; }
	      T& operator()(size_t iii)            { return data_[iii]; }
	const T& v_i(int i, const T* vr)     const { return *(vr+i); }
	      T& v_i(int i, T* vr)                 { return *(vr+i); }
	const T& v_j(int j, const T* vr)     const { return *(vr+j*nnn_.x); }
	      T& v_j(int j, T* vr)                 { return *(vr+j*nnn_.x); }
	const T& v_k(int k, const T* vr)     const { return *(vr+k* nij_); }
	      T& v_k(int k, T* vr)                 { return *(vr+k*nij_); }
	const T* p_i(int i, const T* vr)     const { return (vr+i); }
	const T* p_j(int j, const T* vr)     const { return (vr+j*nnn_.x); }
	const T* p_k(int k, const T* vr)     const { return (vr+k*nij_); }
	size_t   I_i(int i, size_t iii)      const { return (iii+i); }
	size_t   I_j(int j, size_t iii)      const { return (iii+j*nnn_.x); }
	size_t   I_k(int k, size_t iii)      const { return (iii+k*nij_); }
	size_t index(int i, int j, size_t k) const { return k*nij_+j*nnn_.x+i; }

	T*       data()           { return data_.data(); }
	T*       begin()          { return data_.data();  }
	T*       end()            { return &*data_.end();    }
	const T* cbegin()   const { return data_.data(); }
	const T* cend()     const { return &*data_.cend();   }
	const T& back()     const { return data_.back();     }

	const int3& size3() const { return nnn_;    }
	int nx()            const { return nnn_.x;  }
	int ny()            const { return nnn_.y;  }
	int nz()            const { return nnn_.z;  }
	long long  nxy()    const { return nij_;    }
	void getSize(int& n1, int& n2, int& n3) const;

};



//! Base class handling image files with different runtime-time data types (float, uchar, int...)
class voxelImageTBase {
 public:
	virtual ~voxelImageTBase() {}
	virtual void write(std::string fileName) const = 0;
	virtual void printInfo() const {}
	virtual std::unique_ptr<voxelImageTBase> copy() const = 0;
	virtual const int3& size3() const = 0;
	virtual const dbl3& dx() const = 0;
	virtual const dbl3& X0() const = 0;
};


//!  3D image data of different compile-time data types (T = float, uchar, int...)
template <typename T>
class voxelImageT: public voxelImageTBase, public voxelField<T>  {

	dbl3	X0_;   //!< origin
	dbl3	dx_;   //!< voxel size

 public:

	voxelImageT(): X0_(0.,0.,0.), dx_(1,1,1) {}

	voxelImageT(int n1, int n2, int n3, T value) //do not remove, the following constructor will be misused in old codes!
	: voxelField<T>( n1,  n2,  n3,  value),  X0_(0.,0.,0.), dx_(1.,1.,1.) {}

	voxelImageT(int3 n, dbl3 dx=dbl3(1.,1.,1.), dbl3 xmin=dbl3(0.,0.,0.), T value=0)
	: voxelField<T>( n.x,  n.y,  n.z,  value), X0_(xmin),dx_(dx) {}

	voxelImageT(const std::string& fileName, int processKeys=1)
	:	X0_(0.,0.,0.),dx_(1,1,1)  { readFromHeader(fileName, processKeys); }

	void readFromHeader(const std::string& fileName, int processKeys=1);

	bool readAscii(std::string fileName);
	void readRLE(std::string fileName);//! For Avizo format

	std::unique_ptr<voxelImageTBase> copy() const { return std::make_unique<voxelImageT<T>>(*this); }

	void cropD(int3 frm,  int3 to,int emptylyrs=0, T eLyrsValue=1, bool verbose=false) ;

	void writeHeader(std::string fileName) const;
	void writeHeader(std::string fileName, int3 iBgn, int3 iEnd) const { voxelField<T>::writeHeader(fileName,iBgn,iEnd,dx_,X0_); }
	void writeNoHdr(std::string fileName) const;
	void write(std::string fileName) const;

	void erodeLayer(int i);
	void rotate(char direction);
	void PointMedian032(int nAdj0, int nAdj1, T lbl0, T lbl1);
	size_t FaceMedian06(int nAdj0,int nAdj1); // obsolete
	void mode(short nNeist, bool verbose=false); // depricated, use standalone version
	void zeroGrad(int nlyr);

	void AND(const voxelImageT& data2);
	void NOT(const voxelImageT& data2);
	void OR(const voxelImageT& data2);
	void XOR(const voxelImageT& data2);
	void maxEq(const voxelImageT& data2);
	void minEq(const voxelImageT& data2);

	void fillHoles(int maxHoleRadius);

	void shrinkPore();
	void growPore();
	void growLabel(T vl);

   void threshold101(T theresholdMin,T theresholdMax);

	void writeAConnectedPoreVoxel(std::string fileName) const;
	template<typename T2> void resetFrom(const voxelImageT<T2>&Values); // See also: readConvertFromHeader resetFromImageT
	void setFrom(const voxelImageT<T>&Values, int n1, int n2, int n3);
	void growBox(int nlyr);
	void shrinkBox(int nlyr) {  int3 bgn(nlyr,nlyr,nlyr);  cropD(bgn,this->size3()-bgn);  }

	double volFraction(T vv1,T vv2) const;
	void   printInfo() const;

	const dbl3& X0  ()  const { return X0_; }
	dbl3&       X0Ch()        { return X0_; }
	const dbl3& dx  ()  const { return dx_; }
	dbl3&       dxCh()        { return dx_; }
	const int3& size3() const { return voxelField<T>::size3(); }

	double vv_mp5(double i, double j, double k) const { /// linear interpolation
		///set i,j,k -=0.5 (=_mp5) before passing them here, assuming vxl centres are at +0.5
		const int i0=std::min(int(i),this->nx()-2), j0=std::min(int(j),this->ny()-2), k0=std::min(int(k),this->nz()-2); ///. note neg fracs round to zero
		double dd=i-i0,    ld=1.-dd;
		const T*
		vp=&(*this)(i0,j0,k0);  const double v00= *vp*ld + dd*this->v_i(1,vp);
		vp=this->p_j( 1,vp);    const double v10= *vp*ld + dd*this->v_i(1,vp);
		vp=this->p_k( 1,vp);    const double v11= *vp*ld + dd*this->v_i(1,vp);
		vp=this->p_j(-1,vp);    const double v01= *vp*ld + dd*this->v_i(1,vp);
		dd=j-j0;    ld=1.-dd;       k-=k0;
		return (v00*ld + dd*v10) *(1.-k)+k* (v01*ld + dd*v11);
	}
};


template<typename T> void mode26(voxelImageT<T>& vImg, short minDif, bool verbose=false);

std::string VxlKeysHelp(std::string keyname="", std::string subkey="");


template <typename T> using VxlFunc  = bool(*)(std::stringstream&, voxelImageT<T>&);
template <typename T> using VxlFuncs = std::unordered_map<std::string, VxlFunc<T>>;

template<class InpT, typename T>
int vxlProcess(const InpT& inks, voxelImageT<T>& img, std::string nam="");// InpT= string or InputFile

template<class InpT, typename First=uint8_t, typename... Rest>
int vxlProcess(const InpT& inks, voxelImageTBase* ptr, std::string nam="");


std::unique_ptr<voxelImageTBase> readImage(std::string hdrNam /*headername or image type*/, int procesKeys = 1 );

template<class T, typename First=uint8_t, typename... Rest>
int resetFromImageT(voxelImageT<T>& vImg, const voxelImageTBase* imgPtr) { //! cast to specified vImg type
	if(auto img = dynamic_cast<const voxelImageT<First>*>(imgPtr)) { vImg.resetFrom(*img); return 0; }
	else if(sizeof...(Rest)) return resetFromImageT<T,Rest...>(vImg, imgPtr);
	alert("Unknown image type in resetFromImageT");
	return -1;
}

template<typename T>
void readConvertFromHeader( voxelImageT<T>& vImg, std::string hdrNam, int procesKeys=1)  { //! read image and if needed convert its type
	std::unique_ptr<voxelImageTBase> vImgUptr = readImage(hdrNam,procesKeys);
	voxelImageTBase* imgPtr = vImgUptr.get();
	if  (auto img = dynamic_cast<voxelImageT<T>*>(imgPtr)) { vImg = std::move(*img); } //TODO ensure use swap
	else  ensure(  resetFromImageT(vImg, imgPtr)==0, "can not convert image", -1);
}

template<class T> voxelImageT<T>*
vxlCast(void* imgPtr) { return dynamic_cast<voxelImageT<T>*>(static_cast<voxelImageTBase*>(imgPtr)); }


template<typename T>
voxelImageT<T> copyOrReadImgT(std::string hdrNam) {
	#ifdef _STOR_PUB
	if (void* ptr = dbget(_STOR,hdrNam,0)) {
		auto imgPtr = vxlCast<T>(ptr); ensure(imgPtr, "wrong image type", -1);
		return *imgPtr;
	}
	else
	#endif
		return  voxelImageT<T>(hdrNam);
}

#ifdef _STOR_PUB
	#define _dbgetOrReadImg(_T,_img,_hdrNam)  \
			voxelImageT<_T>* _img##Ptr;   std::unique_ptr<voxelImageTBase> _img##Read;  \
			if (void* ptr = dbget(_STOR,_hdrNam,0)) { _img##Ptr = vxlCast<_T>(ptr); }  \
			else       { _img##Read = readImage(_hdrNam,1); _img##Ptr=vxlCast<_T>(_img##Read.get()); }  \
			ensure(_img##Ptr, "wrong image type", -1);  \
			const voxelImageT<_T>& _img = *_img##Ptr;
#else
	#define _dbgetOrReadImg(_T, _img,_hdrNam) voxelImageT<_T> _img(_hdrNam)
#endif

typedef voxelImageT<unsigned char> voxelImage;   //! default image type



#ifndef voxelImageMacros_H
#define voxelImageMacros_H



#define forAlliii_(_vxls)   OMPFor()	\
	for ( size_t iii=0; iii<(_vxls).data_.size(); ++iii)

#define forAllvp_(_vxls)   OMPFor()	\
	for(auto vp=(_vxls).data_.data(); vp<&(*(_vxls).data_.cend()); ++vp)
#define forAll_vp_(_vxls)   OMPFor()	\
	for(auto vp=(_vxls).data_.data()+(_vxls).data_.size()-1; vp>&(*(_vxls).data_.cbegin())-1; --vp)


#define forAllkvp_(_vxls)	OMPFor()	\
	for (int k=0; k<(_vxls).nz(); ++k)   \
	for(auto vp=(_vxls).data_.data()+k*(_vxls).nij_, _ve_=(_vxls).data_.data()+(k+1)*(_vxls).nij_; vp<_ve_; ++vp)
#define forAllk_vp_(_vxls)   OMPFor()	\
	for (int k=0; k<(_vxls).nz(); ++k)   \
	for(auto vp=(_vxls).data_.data()+(k+1)*(_vxls).nij_-1, _ve_=(_vxls).data_.data()+k*(_vxls).nij_-1; vp>_ve_; --vp)

//#define forAllkji(_vxls)   for (int k=0; k<(_vxls).nz(); ++k)  for (int j=0; j<(_vxls).ny(); ++j) for (int i=0; i<(_vxls).nx(); ++i)

#define forkji_be(_ib_,_jb_,_kb_, _ie_,_je_,_ke_)   \
	for (int k=_kb_; k<(_ke_); ++k)   \
	for (int j=_jb_; j<(_je_); ++j)   \
	for (int i=_ib_; i<(_ie_); ++i)
#define forAllkji(_vxls) forkji_be(0,0,0, (_vxls).nx(), (_vxls).ny(), (_vxls).nz())

#define forAllkji_(_vxls)     OMPFor()	forAllkji(_vxls)


#define forAllvp_seq(_vxls)   for(auto vp=(_vxls).data_.data(); vp<&(*(_vxls).data_.cend()); ++vp)
#define forAllcp(_vxls)       for(auto cp=_vxls.data_.cbegin(); cp<_vxls.data_.cend(); ++cp)
#define forAllcp_seq_(_vxls)  for(const auto* cp=_vxls.data_.data(), *_ve_=&(*_vxls.data_.cend()); cp<_ve_; ++cp)
#define forAllvv_seq(_vxls)   for(auto vv : _vxls.data_)
#define forAllvr_seq(_vxls)   for(auto& vr : _vxls.data_)
#define forAlliii_seq(_vxls)  for(size_t iii=0; iii<_vxls.data_.size(); ++iii)

#define forAllkji_m_seq(_nNei,_vxls)   \
 	for (int k=_nNei; k<(_vxls).nz()-_nNei; ++k)   \
 	for (int j=_nNei,_je_=(_vxls).ny()-_nNei; j<_je_; ++j)   \
 	for (int i=_nNei,_ie_=(_vxls).nx()-_nNei; i<_ie_; ++i)
#define forAllkji_1(_vxls)                            forAllkji_m_seq(1,_vxls)
#define forAllkji_1_(_vxls)                  OMPFor() forAllkji_m_seq(1,_vxls)
#define forAllkji_1_sum_(_vxls, _redu)  OMPFor(_redu) forAllkji_m_seq(1,_vxls)
#define forAllkji_m_(_nNei,_vxls)            OMPFor() forAllkji_m_seq(_nNei,_vxls)



#define forAllNei(r1,r2) \
	for (int k_nei_rel=r1; k_nei_rel<=r2; ++k_nei_rel) \
	for (int j_nei_rel=r1; j_nei_rel<=r2; ++j_nei_rel) \
	for (int i_nei_rel=r1; i_nei_rel<=r2; ++i_nei_rel)
#define _nei(_vxls,i_M,j_M,k_M) (_vxls)(i_M+i_nei_rel, j_M+j_nei_rel, k_M+k_nei_rel)

#define forAllNeiInt(ri1,ri2,rj1,rj2,rk1,rk2) \
	for (int k_nei_abs = rk1; k_nei_abs <= rk2; ++k_nei_abs) \
	for (int j_nei_abs = rj1; j_nei_abs <= rj2; ++j_nei_abs) \
	for (int i_nei_abs = ri1; i_nei_abs <= ri2; ++i_nei_abs)
#define _neiv(_vxls) (_vxls)(i_nei_abs, j_nei_abs, k_nei_abs)


//not used in libvoxel:

#define forAlliii_k(_vxls)		\
	for(size_t iii=k*_vxls.nij_; iii<size_t(k+1)*_vxls.nij_; ++iii)

#define forkjid_seq(iMin_m,iMax_m,delta_i)   \
	for (int k=iMin_m.z; k<iMax_m.z; k+=delta_i )   \
	for (int j=iMin_m.y; j<iMax_m.y; j+=delta_i  )   \
	for (int i=iMin_m.x; i<iMax_m.x; i+=delta_i )


#define forAllNeiU(ri1,ri2,rj1,rj2,rk1,rk2) \
	for (int k_nei_m = rk1; k_nei_m <= rk2; ++k_nei_m) \
	for (int j_nei_m = rj1; j_nei_m <= rj2; ++j_nei_m) \
	for (int i_nei_m = ri1; i_nei_m <= ri2; ++i_nei_m)

#define forAllNeiHet(ii1,ii2,jj1,jj2,kk1,kk2) \
	for (short k_nei_m=kk1; k_nei_m<=kk2; ++k_nei_m) \
	for (short j_nei_m=jj1; j_nei_m<=jj2; ++j_nei_m) \
	for (short i_nei_m=ii1; i_nei_m<=ii2; ++i_nei_m)



#define forAllNei_m_(_nNei,r1,r2) \
	for (int k_nei_rel=r1; k_nei_rel<=r2; k_nei_rel+=_nNei) \
	for (int j_nei_rel=r1; j_nei_rel<=r2; j_nei_rel+=_nNei) \
	for (int i_nei_rel=r1; i_nei_rel<=r2; i_nei_rel+=_nNei)

#endif  /// voxelImageMacros_H


// backward compatibility with ICL-2021 versions TODO: delete
#ifndef MS_VERSION
#define MS_VERSION 1
#endif
#if MS_VERSION < 2
#define _ExtraVxlTypes
#include "voxelImageI.h"
#endif
