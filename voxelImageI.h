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



#include "voxelImage.h"
#include <map>
#include <limits>   // std::numeric_limits

#ifdef ZLIB
#include "zfstream.h"
#endif

#ifdef TIFLIB
#include "voxelTiff.h"
#endif


template<typename T> Tint sum6Nei(const voxelImageT<T>& vf, const T* vp)  {
	return Tint(vf.v_i(-1,vp))+vf.v_i( 1,vp)+  vf.v_j(-1,vp)+vf.v_j( 1,vp)+ vf.v_k(-1,vp)+vf.v_k( 1,vp); }

template<typename T> int nSam6Nei(const voxelImageT<T>& vf, const T* vp)  {/// \return number of adjacent voxels with same value, out of 6
	return Tint(vf.v_i(-1,vp)==*vp)+ (vf.v_i( 1,vp)==*vp)+ (vf.v_j(-1,vp)==*vp)+
		(vf.v_j( 1,vp)==*vp)+ (vf.v_k(-1,vp)==*vp)+(vf.v_k( 1,vp)==*vp); }

template<typename T> int nSam6Nei2(const voxelImageT<T>& vf, const T* vp)  {/// \return number of adjacent voxels with same value, out of 6, to add to nSameNei26
	return Tint(vf.v_i(-2,vp)==*vp)+ (vf.v_i( 2,vp)==*vp)+ (vf.v_j(-2,vp)==*vp)+
		(vf.v_j( 2,vp)==*vp)+ (vf.v_k(-2,vp)==*vp)+(vf.v_k( 2,vp)==*vp); }

template<typename T> T maxNei(const voxelImageT<T>& image, int i, int j, int k, int r11, int r22)  {
	T maxx=std::numeric_limits<T>::min();
	forAllNei(r11,r22)  maxx=std::max(maxx,_nei(image,i,j,k));
	return maxx;
}

template<typename T,typename TRes>
TRes accumulate(const voxelField<T>& vf, TRes const & (*operatorFunc)(TRes const &,  TRes const&), TRes result=0)  {
	//result = std::accumulate(vf.data_.begin(), vf.data_.end(), result, operatorFunc);   return result;

	OMPragma("omp parallel")
	{
		TRes res_private = result;
		OMPragma("omp for nowait")
		for(auto vp=vf.data_.begin(); vp<vf.data_.end(); ++vp)
		{
			res_private=operatorFunc(res_private,*vp);
		}
		OMPragma("omp critical")
		result = operatorFunc(res_private,result);
	}
	return result;
}

template<typename T,typename TRes, typename TFunc>
TRes accumulate(const voxelField<T>& vf, TFunc operatorFunc, TRes result=0)  {
	//result = std::accumulate(vf.data_.begin(), vf.data_.end(), result, operatorFunc);   return result;

	OMPragma("omp parallel")
	{
		TRes res_private = result;
		OMPragma("omp for nowait")
		for(auto vp=vf.data_.begin(); vp<vf.data_.end(); ++vp)  {
			res_private=operatorFunc(res_private,*vp);
		}
		OMPragma("omp critical")
		result = operatorFunc(res_private,result);
	}
	return result;
}

template<typename T>
int accumulateT(const voxelImage& vf, std::function<T(T, T)> operatorFunc, T result=0)  {
	return std::accumulate(vf.data_.begin(), vf.data_.end(), result, operatorFunc);
}



#define dblfunc(M_func_) (double const & (*) (double const &, double const &))(M_func_<double>)


template<class voxelFieldT>
vars<dbls> vxlDist(const voxelFieldT& vf, int nsteps=32, double minV=3e38, double maxV=-3e38)
{
	vars<dbls> distrib(2,dbls(nsteps,0.));

	if(minV> 1e38) minV=accumulate(vf,dblfunc(std::min), 1e38);
	if(maxV<-1e38) maxV=accumulate(vf,dblfunc(std::max),-1e38);
	double spanV=(maxV-minV);
	double deltaV=spanV/nsteps*(1.+1e-14);
	std::cout<<"  vxlDist_range: ["<<minV<<" "<<maxV<<"]";

	for (int i=0; i<nsteps; ++i)	distrib[0][i] = minV+deltaV/2+i*deltaV;

	for (auto vv : vf.data_ ) distrib[1][std::min(std::max(0, int((vv-minV)/deltaV+0.5)),nsteps)]+=1.;
	ensure(distrib[1].sum()>1e-32);
	distrib[1]/=distrib[1].sum()*deltaV;

	return distrib;
}




template<typename T>   void voxelField<T>::reset(int3 nnn)  {
	nij_=size_t(nnn.x)*nnn.y;
	this->data_.resize(nnn.z*nij_);
	nnn_=nnn;
}

template<typename T>   void voxelField<T>::reset(int3 nnn, T value)  {
	nij_=size_t(nnn.x)*nnn.y;
	this->data_.assign(size_t(nnn.z)*nij_, value);
	nnn_=nnn;
}


template<typename T>   void voxelField<T>::getSize(int& n1, int& n2, int& n3) const {
	n3 = (*this).nz();
	if (n3>0) { n1 = (*this).nx();  n2 = (*this).ny();  }
	else      { n1 = 0;  n2 = 0; }
}






//read order sensitive
template<typename T>   void voxelField<T>::readAscii(std::ifstream& in)  {
	forAllvr_seq((*this)) in>>vr;
}

template<>  inline  void voxelField<unsigned char>::readAscii(std::ifstream& in)  {
	int tmp;
	forAllvr_seq((*this)) { in>>tmp;  vr=tmp; }
}

//read order sensitive
template<typename T>   bool voxelField<T>::readAscii(std::string fnam)  {
	std::cout<<" ascii, reading "<<fnam<<std::endl;
	std::ifstream in(fnam);  ensure(in, "cannot open image file "+fnam,-1);
	readAscii(in);
	return !in.fail();
}

template<typename T>   bool voxelImageT<T>::readAscii(std::string fnam)  {
		//!  overwrite voxelField<T>::readAscii(), as voxelField interprets  8bit numerical values as characters not integers

	std::cout<<" reading "<<fnam<<std::endl;

	std::ifstream in(fnam);
	assert(in);

	char tmpc[8];
	for (int i=0; i<8; ++i)   in>>tmpc[i];
	if (std::string(tmpc).compare(0,4,"ascii") == 0) //ignore first lines
	{
		int3 nnn;           in>>nnn.z>>nnn.x>>nnn.y;//ignore first lines
		dbl3  xmin,xmax;  in>> xmin.x>>xmax.x>>xmin.y>>xmax.y>>xmin.z>>xmax.z ;
		std::cout<<"Warning: ignoring the header of file "<<fnam<<std::endl;
	}
	else
		in.seekg(0, in.beg);
	voxelField<T>::readAscii(in);
	return !in.fail();
}



//read order sensitive
template<typename T>   void voxelField<T>::readMicroCT(std::string fnam)  {
	std::cout<<" micro-CT, reading "<<fnam<<std::endl;
	std::ifstream in(fnam);    ensure(in,"can not open image file, "+fnam,-1);

	char tmpc;
	for (int i=0; i<8; ++i)   in>>tmpc, std::cout<<" "<<tmpc;  //ignore the first 8 characters (ascii 3uc)
	int3 nnn;
	dbl3  xmin,  xmax;

	in>>nnn.z>>nnn.y>>nnn.x;						// number of variables (dimension of
	ensure(nnn.x==this->nx());
	in>>	xmin  >>	xmax ;
	readAscii(in); // no reset?!!
}


inline std::string getAmiraDataType(const std::string& fnam)  {
	std::ifstream hdr(fnam);
	ensure(hdr, "could not open "+fnam, -1);
	std::string tmp;
	while (true)  {
		hdr>>tmp;
		std::stringstream ss;
		if(hdr.peek()!='\n') hdr.get (*(ss.rdbuf()));
		if (hdr.fail()) {std::cout<<"Error reading "<<fnam<<",  after "<<tmp<<std::endl; break;}

		if (tmp=="Content")  { ss >> tmp >> tmp;  break;  }
		else if (tmp=="@1")    break;
	}
	if(tmp.size() && tmp.back()==',') tmp.resize(tmp.size()-1);
	return tmp;
}

inline void getAmiraHeaderSize(const std::string& fnam, int3& nnn, dbl3& dx_, dbl3& X0_, int& nSkipBytes, int& RLE)  {
	std::ifstream hdr(fnam);
	while (true)  {
		std::string tmpStr;    hdr>>tmpStr;

		std::stringstream ss;
		if(hdr.peek()!='\n') hdr.get (*(ss.rdbuf()));
		if (hdr.fail()) {std::cout<<"Error reading "<<fnam<<",  after "<<tmpStr<<std::endl; break;}
		std::string tmp;
		if (tmpStr == "define")  {
			ss >> tmp  >>nnn;
			if (tmp != "Lattice") std::cout<<" Warning: define != Lattice n3, read: "<<tmp<<std::endl;
		}
		else if (tmpStr == "BoundingBox")  {
			ss >> X0_.x>>dx_.x>> X0_.y>>dx_.y>> X0_.z>>dx_.z;
			int3 nn=nnn; for(int i:{0,1,2}) nn[i]=std::max(nn[i]-1,1);// Why nnn-1? well it seems Avizo cannot even properly convert voxel size to bounding box, or something else is wrong in Abdollah way of doing things!!!
			dx_=(dx_-X0_)/nn;
		}
		else if (tmpStr=="@1")     break;
		else if(tmpStr=="Lattice")  {
			while (tmpStr[0]!='@' && ss) 	 ss>>tmpStr;
			RLE = tmpStr.size()>11 && tmpStr.compare(3,9,"HxByteRLE") == 0;
		}
	}
	if (nnn.z>1.1*maxNz+2)  nnn.z = maxNz;
	nSkipBytes = hdr.tellg(); ++nSkipBytes; //++ is for '\n' after "@1"
}

template<typename T>   int voxelField<T>::readBin(std::string fnam, int nSkipBytes)  {
	int3 nnn = size3();
	int RLEcompressed=0;

	#ifdef TIFLIB
	if(hasExt(fnam,".tif"))      return readTif(*this, fnam);
	#endif

	(std::cout<<"  Reading "<<fnam<<" ").flush();

	if(hasExt(fnam,".am"))  {
		dbl3	X0, dx;
		getAmiraHeaderSize(fnam, nnn,dx,X0,nSkipBytes,RLEcompressed);
		(std::cout<<", .am  format").flush();
		this->reset(nnn);
	}

	(std::cout<<", size:"<<size_t(nnn.x)<<'*'<<nnn.y<<'*'<<nnn.z<<"*"<<sizeof(T)).flush();

	if(hasExt(fnam,".gz"))  {

	 #ifdef ZLIB
		if(std::ifstream(fnam).good())  {
			(std::cout<<", using libz").flush();
			gzifstream  in(fnam.c_str());
			in.read(reinterpret_cast<char*>(&((*this)(0,0,0))), (size_t(nnn.x)*nnn.y*nnn.z)*sizeof(T) );
			in.close();

			std::cout<<"."<<std::endl;
			return 0;
		} else std::cout<<"Error: cannot read "<<fnam<<std::endl;
	 #endif

		fnam=fnam.substr(0,fnam.size()-3);
		std::cout<<" .gz not read or not supported, trying "<<fnam<<" instead"<<std::endl;

	}

	std::ifstream in (fnam, std::ios::in | std::ios::binary);
	ensure(in, "cannot open image file "+fnam,-1);
	if(nSkipBytes) in.ignore(nSkipBytes);

	if(RLEcompressed)  {
		std::cout<<", RLE decoding";
		char count=0;
		char val=0;
		//ensure(sizeof(T)==1,"Only 8bit .am files are supported");
		char*      vp = reinterpret_cast<char*>(voxelField<T>::data_.data());
		char*const ve =reinterpret_cast<char*>(std::to_address(voxelField<T>::data_.end()));
		while(vp<ve)  {
			in.get(count); in.get(val);
			if(count&char(0x80)) {
				*vp=val;
				count&=0x7f;
				if (count > ve-vp) [[unlikely]] count = char(ve-vp); // for maxNz
				while((--count)!=0) { in.get(val); *(++vp)=val; }
				++vp;
			}
			else  {
				if (count > ve-vp) [[unlikely]] count = char(ve-vp); // for maxNz
				std::fill(vp,vp+count,val);  vp+=count;
			}
		}
	}
	else  {
		(std::cout<<  ", reading raw data").flush();
		//in.read(reinterpret_cast<char*>(&((*this)(0,0,0))), (size_t(nnn.x)*nnn.y*nnn.z)*sizeof(T) ); // fails on crappy MSWindows
		for (int k=0; k<nnn.z;++k)
			in.read(reinterpret_cast<char*>(&((*this)(0,0,k))), (size_t(nnn.x)*nnn.y)*sizeof(T));
	}
	(std::cout<<". ").flush();

	if (!in)  {  std::cout<<"\n\n ***** Error in reading "<<fnam<<" ***** \n"<<std::endl;  return -1;  }
	else return 0;

}



template<typename T>
int voxelField<T>::readBin(std::string fnam, int iS,int iE, int jS,int jE, int kS,int kE, int nSkipBytes)  {
	if(hasExt(fnam,".tif"))  {
		voxelImageT<T> vxls(fnam);
		if (vxls.nx()!=iE-iS)	std::cout<<"Error in reading "<<fnam<<", unexpected size: Nx="<<vxls.nx()<<"  !="<<iE-iS<<std::endl;
		setBlock(iS,jS,kS, vxls);
		return 0;
	}
	if(hasExt(fnam,".gz"))  {
		voxelField<T> vxls(int3(iE-iS, jE-jS, kE-kS));
		vxls.readBin(fnam);
		setBlock(iS,jS,kS, vxls);
		return 0;
	}

	(std::cout<<" reading binary file "<<fnam<<" ").flush();
	std::ifstream in (fnam, std::ios::in | std::ios::binary);
	ensure(in, "cannot open image file "+fnam,-1);
	if(nSkipBytes) in.ignore(nSkipBytes);
	int k = kS;
	for ( ; k < kE; k++)  {
		for (int j = jS; j<jE; ++j)
			if(in)  in.read(reinterpret_cast<char*>(&((*this)(iS,j,k))), (iE-iS)*sizeof(T) );
		if (!in) break;
	}
	if (!in)	std::cout<<"\n\n ***** Error in reading "<<fnam<<" ***** \n"<<"only "<<k<<" layers read"<<std::endl;
	std::cout<<std::endl;
	return 0;
}


template<typename T>   void voxelField<T>::writeBin(std::string fnam) const  {
	int3 nnn = size3();

	if(hasExt(fnam,".tif"))  {
	 #ifdef TIFLIB
		(std::cout<<"  writing tif file "<<fnam<<";  size: "<<nnn<<" ").flush();
		writeTif(*this, fnam);
		std::cout<<std::endl;
		return;
	 #else
		fnam = fnam.substr(0,fnam.size()-4)+imgExt();
	 #endif //TIFLIB
	}

	#ifdef ZLIB
	if(hasExt(fnam,".gz"))  {
		(std::cout<<"  writing compressed file "<<fnam<<";  size: "<<nnn).flush();
		gzofstream  of(fnam.c_str());
		of << setcompression(5);//,Z_RLE TODO: benchmark  Z_DEFAULT_COMPRESSION   Z_RLE +gz
		ensure(of,"Failed in writing "+fnam);
		if(data_.size())
		of.write(reinterpret_cast<const char*>(&((*this)(0,0,0))), (size_t(nnn.x)*nnn.y*nnn.z) * sizeof(T));
		of.flush();
		std::cout<<std::endl;
		return;
	}
	#else
	if(hasExt(fnam,".gz")) fnam=fnam.substr(0,fnam.size()-3);
	#endif


	(std::cout<<"  writing binary file "<<fnam<<";  size: "<<nnn).flush();
	auto mod = std::ios::out | std::ios::binary;
	if( hasExt(fnam,".am")) 	{ // in case write Header called before
		char cs[]="xxx";
		std::ifstream is(fnam); if(is)	{ is.seekg (3, is.end);	is.get(cs,3); }		is.close();
		if(cs[0]!='@' || cs[1]!='1' || cs[2]!='\n')	writeHeader(fnam,{0,0,0},nnn,{1.,1.,1.},{0.,0.,0.});
		mod = mod | std::ios::app;
	}
	//else if(!hasExt(fnam,".tif")) writeHeader(fnam,{0,0,0},nnn); this overwrites voxelImage one

	std::ofstream of (fnam, mod);
	assert(of);
	if(data_.size())
	//of.write(reinterpret_cast<const char*>(&((*this)(0,0,0))), (size_t(nnn.x)*nnn.y*nnn.z) * sizeof(T)); // fails on crappy MSWindows
	for (int k=0; k<nnn.z;++k)
		of.write(reinterpret_cast<const char*>(&((*this)(0,0,k))), (size_t(nnn.x)*nnn.y)*sizeof(T));
	of.flush();
	std::cout<<std::endl;

}


template<typename T>
void voxelField<T>::writeBin(std::string fnam, int iS,int iE, int jS,int jE, int kS,int kE ) const  {


	if(hasExt(fnam,".tif"))  {
	 #ifdef TIFLIB
		(std::cout<<"  writing tif file "<<fnam<<";  i: "<<iS<<" "<<iE<<",  j: "<<jS<<" "<<jE<<",  k: "<<kS<<" "<<kE).flush();
		writeTif(*this, fnam, iS, iE,  jS, jE,  kS, kE);
		std::cout<<std::endl;
		return;
	 #else
		fnam = fnam.substr(0,fnam.size()-4)+imgExt();
	 #endif //TIFLIB
	}

	#ifdef ZLIB
	if(hasExt(fnam,".gz"))  {
		(std::cout<<"  writing compressed file "<<fnam<<";  i: "<<iS<<" "<<iE<<",  j: "<<jS<<" "<<jE<<",  k: "<<kS<<" "<<kE).flush();
		gzofstream  of(fnam.c_str());
		of << setcompression(5);//,Z_RLE
		assert(of);
		if(data_.size())
		 for (int k = kS; k < kE; k++)
			for (int j = jS; j < jE; ++j)
			{
				of.write(reinterpret_cast<const char*>(&((*this)(iS,j,k))), (iE-iS) * sizeof(T));
			}
		of.flush();
		of.close();
		std::cout<<std::endl;
		return;
	}
	#endif

	(std::cout<<"  writing binary file "<<fnam<<";  i: "<<iS<<" "<<iE<<",  j: "<<jS<<" "<<jE<<",  k: "<<kS<<" "<<kE).flush();
	auto mod = std::ios::out | std::ios::binary;
	if(hasExt(fnam,".am"))  {
		char cs[]="xxx";
		std::ifstream is(fnam);		if(is)	{ is.seekg (3, is.end);	is.get(cs,3); }		is.close();
		if(cs[0]!='@' || cs[1]!='1' || cs[2]!='\n')	writeHeader(fnam,{iS,jS,kS},{iE,jE,kE}, {1.,1.,1.}, {0.,0.,0.});
		mod = mod | std::ios::app;
	}
	//else if(!hasExt(fnam,".tif"))		writeHeader(fnam,{iS,jS,kS},{iE,jE,kE}, 1., 0.); this overwrites voxelImage one
	std::ofstream of(fnam, mod);
	assert(of);
	if(data_.size())
	 for (int k = kS; k < kE; k++)
		for (int j = jS; j < jE; ++j)
			of.write(reinterpret_cast<const char*>(&((*this)(iS,j,k))), (iE-iS) * sizeof(T));
	std::cout<<std::endl;
}

template<typename T>
void voxelField<T>::writeAscii(std::string fnam,int iS,int iE, int jS,int jE, int kS,int kE) const  {
	(std::cout<<"  writing ascii file "<<fnam<<";  ").flush();

	std::ofstream of (fnam);
	assert(of);

	for (int k = kS; k < kE; k++)
		for (int j = jS; j < jE; ++j) {
			for (int i = iS; i < iE; ++i)
				of<<((*this)(i,j,k))<<' ';
			of<<"\n";
		}
	of<<std::endl;
	std::cout<<std::endl;
}

template<> inline
void voxelField<unsigned char>::writeAscii(std::string fnam,int iS,int iE, int jS,int jE, int kS,int kE) const  {
	(std::cout<<"  writing ascii file "<<fnam<<";  ").flush();

	std::ofstream of (fnam);
	assert(of);

	for (int k = kS; k < kE; k++)
	  for (int j = jS; j < jE; ++j)
	  {
		for (int i = iS; i < iE; ++i)
		   of<<int((*this)(i,j,k))<<' ';
		of<<"\n";
	  }
	of<<std::endl;
	std::cout<<std::endl;
}

template<typename T>   void voxelField<T>::writeAscii(std::string fnam) const  {
	int3 imgsize = size3();
	writeAscii(fnam,0, imgsize.x, 0,imgsize.y, 0,imgsize.z);
}

template<typename T>   void voxelField<T>::writeRotatedXZ(std::ofstream& of) const  {
	for (int i=0; i<(*this).nx(); ++i) //reversed order with k
		for (int j=0; j<(*this).ny(); ++j)  {
			for (int k=0; k<(*this).nz(); k++) //reversed order with i
				of<<double((*this)(i,j,k))<<' ';
			of<<std::endl;
		}
	std::cout<<"  write RotatedXZ. "<<std::endl;
}



template<typename T> void voxelImageT<T>::writeHeader(std::string ofnam) const  {
	voxelField<T>::writeHeader(ofnam, int3(0,0,0), voxelField<T>::size3(),dx_,X0_);
}

template<typename T>
void voxelField<T>::writeHeader(std::string fnam, int3 iS, int3 iE, dbl3 dx, dbl3 X0) const  {

	if(dx.x<0)  {
		std::cerr<<"Error negative dx, writing abs value instead" ;
		dx.x=std::abs(dx.x);  dx.y=std::abs(dx.y);  dx.z=std::abs(dx.z);
	}
	dbl3 xmin(X0+iS*dx);
	dbl3 xmax(X0+iE*dx);
	int3	nnn(iE-iS);
	if (hasExt(fnam,".am"))  {
		std::string typ="uchar";
		if (typeid(T)==typeid(char)) typ="char";
		else if (typeid(T)==typeid(short)) typ="short";
		else if (typeid(T)==typeid(unsigned short)) typ="ushort";
		else if (typeid(T)==typeid(int)) typ="int";
		else if (typeid(T)==typeid(int)) typ="uint";
		else if (typeid(T)==typeid(float)) typ="float";
		else if (typeid(T)==typeid(double)) typ="double";
		else if (typeid(T)==typeid(float3)) typ="float[3]";

		std::ofstream out(fnam);
		assert(out);
		out <<"# Avizo BINARY-LITTLE-ENDIAN 2.1\n\n\n"
		    <<"define Lattice "<<nnn<<"\n\n"
		    <<"Parameters {\n    Units {\n        Coordinates \"m\"\n    }\n";
			if(typ=="float[3]")
				out <<"    XLabExperiment {\n        viscosity 0.001,\n        inputPressure 1,\n        outputPressure 0,\n        flowRate 1\n    }\n";

		out <<"    Content \""<<nnn.x<<"x"<<nnn.y<<"x"<<nnn.z<<" "<<typ<<", uniform coordinates\",\n"
		    <<"    BoundingBox "<<xmin.x<<" "<<xmax.x<<" "<<xmin.y<<" "<<xmax.y<<" "<<xmin.z<<" "<<xmax.z<<",\n"
		    <<"    CoordType \"uniform\"\n}\n\n"
		    <<"Lattice { "<<typ<<" Data } @1\n\n# Data section follows\n@1\n";

	} else
	if (hasExt(fnam,"_header"))  {
		std::ofstream of(fnam); assert(of);
		of
		 <<"Nxyz\n"
		   "dxX0\n"
		 <<nnn<<'\n'<<dx<<'\n'<<xmin<<std::endl
		 <<"\n\nComments:\n"
		   " first 9 entries above are:"
		   "	Nx Ny Nz\n	dx dy dz\n	Xo Yo Zo\n"
		   " Nx, Ny and Nz  count for the number of columns, rows and layers respectively as written in the file\n"
		   " Optional keywords (move above Comments to activate):\n"
		   "	crop		0  299   0  299   0  299\n"
		   "	pore 		0 0\n"
		   "	resample	1\n"
		   "	direction	z\n"
		   "	..... "
		 <<std::endl;
	}
	else
	{
		int islash=fnam.find_last_of("\\/"); if (islash>=int(fnam.size())) islash=-1;
		std::string title=fnam.substr(islash+1);
		if (hasExt(fnam,".mhd"))    title=title.substr(0,title.size()-4)+imgExt();
		else if (hasExt(fnam,".raw.gz"))  fnam=fnam.substr(0,fnam.size()-7)+".mhd";
		else  fnam=fnam.substr(0,fnam.find_last_of("."))+".mhd";


		std::string typeNmeVTK="MET_UCHAR";
		if (typeid(T)==typeid(char)) typeNmeVTK="MET_CHAR";
		else if (typeid(T)==typeid(short)) typeNmeVTK="MET_SHORT";
		else if (typeid(T)==typeid(unsigned short)) typeNmeVTK="MET_USHORT";
		else if (typeid(T)==typeid(int)) typeNmeVTK="MET_INT";
		else if (typeid(T)==typeid(int)) typeNmeVTK="MET_UINT";
		else if (typeid(T)==typeid(float)) typeNmeVTK="MET_FLOAT";
		else if (typeid(T)==typeid(double)) typeNmeVTK="MET_DOUBLE";

		std::ofstream of(fnam);  ensure(of, fnam);
		of
			 <<"ObjectType =  Image"<<std::endl
			 <<"NDims =	   3"<<std::endl
			 <<"ElementType = "<<typeNmeVTK <<std::endl
			 <<"ElementByteOrderMSB = False\n"//(typeNmeVTK=="MET_UCHAR" ? "\n":)
			 <<"ElementNumberOfChannels = 1\n"//(typeNmeVTK=="MET_UCHAR" ? "\n":)
			 <<"CompressedData = "<<(hasExt(title,".gz") ? "True" : "False")<<std::endl
			 <<"\nDimSize =   "<<nnn<<std::endl
			 <<"ElementSize =	"<<dx<<std::endl
			 <<"Offset =      "<<X0<<std::endl
			 <<"ElementDataFile = "<<title<<std::endl <<std::endl;
			 if(dx.x>=0.001)
				of <<"Unit = "<<1<<std::endl;

			of <<std::endl <<std::endl;
	}
}


template<typename T> void voxelField<T>::writeNoHdr(std::string fnam) const  {
	if (hasExt(fnam,".mhd")) fnam = fnam.substr(0,fnam.size()-4)+imgExt();
	if (hasExt(fnam,".dat") || hasExt(fnam,".txt"))  this->writeAscii(fnam);
	else if(fnam!="NO_WRITE") this->writeBin(fnam);
}

template<typename T> void voxelImageT<T>::writeNoHdr(std::string fnam) const  {
	if(hasExt(fnam,".am")) writeHeader(fnam);// Warning: remember to append to this in writeBin
	voxelField<T>::writeBin(fnam);
}

template<typename T> void voxelImageT<T>::write(std::string fnam) const  {

	if (hasExt(fnam,".dat") || hasExt(fnam,".txt"))  {
		this->writeAscii(fnam);
		writeHeader(fnam+"_header");
	}
	else if(fnam!="NO_WRITE")  {
		if (hasExt(fnam,".mhd"))  writeHeader(fnam = fnam.substr(0,fnam.size()-4)+imgExt());
		else if(!hasExt(fnam,".tif"))	writeHeader(fnam);//&& !hasExt(fnam,".am")!=0

		this->writeBin(fnam);
	}
}










template<typename T>
void vxlCopy(const voxelImageT<T>& img2, int3 frm,  int3 to, voxelImageT<T>& img, int3 at)  {
	for (int k=frm.z; k<to.z; k++)
		for (int j=frm.y; j<to.y; ++j)
			std::copy(&img2(frm.x, j, k),  &img2(to.x, j, k),  &img(at.x, at.y+j-frm.y, at.z+k-frm.z));
}





template<typename T>
void voxelImageT<T>::cropD( int3 frm,  int3 to, int emptylyrs, T eLyrsValue, bool verbose)  {
	//crop(cropBgn.x,cropEnd.x-1,cropBgn.y,cropEnd.y-1,cropBgn.z,cropEnd.z-1,emptylyrs,eLyrsValue);
	{
		if(to.x<=0) to.x = this->nx()+to.x;
		if(to.y<=0) to.y = this->ny()+to.y;
		if(to.z<=0) to.z = this->nz()+to.z;
	}

	if(verbose) (std::cout<<"  cropping from  ["<<frm <<" to "<<to<<  ")  " ).flush();
	else (std::cout<<" D " ).flush();
	ensure(to.x<=size3().x && to.y<=size3().y && to.z<=size3().z,"cropping outside bounds: "+_s(size3()),2);

	X0_+=(frm-emptylyrs)*dx_;

	voxelImageT<T> tmp;
	tmp.swapData(*this);

	if (emptylyrs)
	{
		if(verbose) (std::cout<<", adding "<<emptylyrs<<" layers of "<< double(eLyrsValue)<<"  ").flush();
		this->reset(to-frm+2*emptylyrs, eLyrsValue);
	}
	else this->reset(to-frm);


	vxlCopy(tmp, frm, to, *this, int3(emptylyrs,emptylyrs,emptylyrs));
}




template<typename T>  void voxelField<T>::setSlice(char dir, int ijk, T vv)  {
	if(dir=='i')
	 for ( size_t iii=ijk; iii<(*this).data_.size() ; iii+=(*this).nx())
		(*this).data_[iii]=vv;
	else if(dir=='j')
	 for (int k=0; k<(*this).nz(); ++k)
		std::fill( &(*this)(0,ijk,k), &(*this)(0,ijk+1,k),vv);
	else if(dir=='k')
		std::fill( &(*this)(0,0,ijk), &(*this)(0,0,ijk+1),vv);
	else 	std::cout<<"Error: wrong dir "<<dir<<std::endl;
}

template<typename T>  void voxelField<T>::setLayer(int k, const T* Values)  {
	//(*this)[k]=Values;
	std::copy(Values, Values+(*this).nx()*(*this).ny(), &(*this)(0,0,k));

}
template<typename T>  void voxelField<T>::replacezLayer(int k, int fromj)  {
	this->setLayer(k,&(*this)(0,0,fromj));
}

template<typename T>  void voxelField<T>::replacexLayer(int i, int fromi)  {

	 for (int k=0; k<(*this).nz(); ++k)
		for (int j=0; j<(*this).ny(); ++j)
				(*this)(i,j,k)=(*this)(fromi,j,k);

}
template<typename T>  void voxelField<T>::replaceyLayer(int j, int fromj)  {

	for (int k=0; k<(*this).nz(); ++k)
			for (int i=0; i<(*this).nx(); ++i)
				(*this)(i,j,k)=(*this)(i,fromj,k);
}


template<typename T>
void voxelField<T>::setBlock(int n1, int n2, int n3, const voxelField<T>& Values)  {
	forAllkji_(Values)
			(*this)(i+n1,j+n2,k+n3)=Values(i,j,k);
}
template<typename T>
void voxelField<T>::setFrom(const voxelField<T>&Values, int n1, int n2, int n3)  { // from image with  lager size,
	forAllkji_(*this)	(*this)(i,j,k)=Values(i+n1,j+n2,k+n3);
}

template<typename T>
void voxelField<T>::swapData(voxelField<T>&other) noexcept {
	std::swap(nij_, other.nij_);
	std::swap(nnn_, other.nnn_);
	std::swap(data_, other.data_);
}

template<typename T>
template<typename T2>
void voxelImageT<T>::resetFrom(const voxelImageT<T2>&Values)  {
	dx_= Values.dx();
	X0_ = Values.X0();
	this->reset(Values.size3(),T(0));
	forAlliii_((*this))
		(*this)(iii)=T(Values(iii));
}
template<typename T>
void voxelImageT<T>::setFrom(const voxelImageT<T>&Values, int n1, int n2, int n3)  { // from image with  lager size
	dx_= Values.dx();
	X0_.x = Values.X0().x+n1*dx_.x;
	X0_.y = Values.X0().y+n2*dx_.y;
	X0_.z = Values.X0().z+n3*dx_.z;
	forAllkji_(*this)
			(*this)(i,j,k)=Values(i+n1,j+n2,k+n3);
}

// obselete, use growBounds instead
template<typename T> void voxelImageT<T>::growBox(int nLayers) {

	int3 nnn = (*this).size3();
	(*this).cropD(int3(0,0,0),nnn, nLayers,1,false);

	for (int i=0; i<nLayers; ++i)  {
		(*this).replaceyLayer(nnn.y+nLayers+i, nnn.y+nLayers-1);
		(*this).replaceyLayer(i, nLayers);
		(*this).replacexLayer(nnn.x+nLayers+i, nnn.x+nLayers-1);
		(*this).replacexLayer(i, nLayers);
		(*this).setLayer(nnn.z+nLayers+i, &(*this)(0,0,(nnn.z+nLayers-1)));
		(*this).setLayer(i, &(*this)(0,0,nLayers));
	}
}


template<typename T>
void voxelImageT<T>::zeroGrad(int nLayers)  {

	int3 nnn = (*this).size3();
	for (int i=0; i<nLayers; ++i)  {
		(*this).replaceyLayer(nnn.y-nLayers+i, nnn.y-nLayers-1);
		(*this).replaceyLayer(i, nLayers);
		(*this).replacexLayer(nnn.x-nLayers+i, nnn.x-nLayers-1);
		(*this).replacexLayer(i, nLayers);
		(*this).setLayer(nnn.z-nLayers+i, &(*this)(0,0,(nnn.z-nLayers-1)));
		(*this).setLayer(i, &(*this)(0,0,nLayers));
	}
}

template<typename T>
voxelImageT<T> growBounds(const voxelImageT<T>& vxls, int nLayers)  {
	int3 nnn = vxls.size3();

	voxelImageT<T> tmp;

	tmp.reset(nnn+2*nLayers);

	for (int k=0; k<nnn.z; k++)
		for (int j=0; j<nnn.y; ++j)
			std::copy(&vxls(0,j+0,k+0), &vxls(0,j+0,k+0)+nnn.x, &tmp(nLayers,j+nLayers,k+nLayers));


	for (int i=0; i<nLayers; ++i)  {
		tmp.replaceyLayer(nnn.y+nLayers+i, nnn.y+nLayers-1);
		tmp.replaceyLayer(i, nLayers);
		tmp.replacexLayer(nnn.x+nLayers+i, nnn.x+nLayers-1);
		tmp.replacexLayer(i, nLayers);
		tmp.setLayer(nnn.z+nLayers+i, &tmp(0,0,(nnn.z+nLayers-1)));
		tmp.setLayer(i, &tmp(0,0,nLayers));
	}
	return tmp;
}


template<typename T>  void voxelImageT<T>::growLabel(T vl)  {

	const voxelImageT<T> vxls=growBounds(*this,1);

	forAllkji_1_(vxls)  {
		if (vxls(i,j,k)==vl)  {
			const T* optr=&vxls(i,j,k);
			      T* vptr=&(*this)(i-1,j-1,k-1);
			if(vxls.v_i(-1,optr)!=vl) (*this).v_i(-1,vptr)=vl;
			if(vxls.v_i( 1,optr)!=vl) (*this).v_i( 1,vptr)=vl;
			if(vxls.v_j(-1,optr)!=vl) (*this).v_j(-1,vptr)=vl;
			if(vxls.v_j( 1,optr)!=vl) (*this).v_j( 1,vptr)=vl;
			if(vxls.v_k(-1,optr)!=vl) (*this).v_k(-1,vptr)=vl;
			if(vxls.v_k( 1,optr)!=vl) (*this).v_k( 1,vptr)=vl;
		}
	}
}


template<typename T>
voxelImageT<T>  resampleMean(const voxelImageT<T>& img, double nReSampleNotSafe) {//  TODO to be tested

	voxelImageT<T> tmp;
	if (nReSampleNotSafe < .999)  {
		int nReS=1./nReSampleNotSafe+0.5;
		tmp.reset(nReS*img.size3());
		forAllkji_(tmp)
			tmp(i,j,k)=img((0.5+i)/nReS,(0.5+j)/nReS,(0.5+k)/nReS);

		tmp.dxCh()=img.dx()/nReS; //fixed
		tmp.X0Ch()=img.X0()/nReS; //fixed
		return tmp;
	}
	else if (nReSampleNotSafe > 1.001)  {
		int nReS=nReSampleNotSafe+0.5; /// Warning unsigned doesn't work  wTf
		tmp.reset(img.size3()/nReS);
		forAllkji_(tmp)  {
			int neiSum=0;
			forAllNei(0,nReS-1) 	neiSum+=_nei(img,i*nReS,j*nReS,k*nReS);
			tmp(i,j,k)=(0.5+double(neiSum)/(nReS*nReS*nReS));
		}
		tmp.dxCh()=img.dx()*nReS; //fixed
		tmp.X0Ch()=img.X0()*nReS; //fixed
		return tmp;
	}
	else return img;
}


template<typename T>
class mapComparer  {  public: bool operator() (std::pair<const T,short>& i1, std::pair<const T,short> i2) {return i1.second<i2.second;}  };


template<typename T>
voxelImageT<T>  resliceZ(const voxelImageT<T>& img, double nReSampleNotSafe)//  TODO to be tested
{
	if (nReSampleNotSafe<1e-64) { nReSampleNotSafe=img.dx().x/img.dx().z; (std::cout << " nResample: "<<nReSampleNotSafe<<", ").flush(); }
	voxelImageT<T> tmp;
	int3 nnn=img.size3();
	if (nReSampleNotSafe < .999)  {
		//int nReS=1./nReSampleNotSafe+0.5;
		double nReS=nReSampleNotSafe;
		tmp.reset({nnn.x, nnn.y, int(nnn.z/nReS)});
		forAllkji_(tmp)
			tmp(i,j,k)=img(i,j,(0.5+k)*nReS);

		tmp.dxCh()=img.dx(); tmp.dxCh().z*=nReS; //fixed
		tmp.X0Ch()=img.X0(); tmp.X0Ch().z*=nReS; //fixed
		return tmp;
	}
	else if (nReSampleNotSafe > 1.001)  {
		std::cout<<"resliceZ >1 not tested" <<std::endl;
		int nReS=nReSampleNotSafe+0.5; /// Warning unsigned doesn't work  wTf
		tmp.reset({nnn.x, nnn.y, nnn.z/nReS});
		forAllkji_(tmp)  {
			std::map<T,short> neis;///.  ID-counter
			Tint sumv=0;
			for_0(nReS,kk) sumv+=img(i,j,k*nReS+kk);
			tmp(i,j,k)=sumv/nReS;
		}
		tmp.dxCh()=img.dx()*nReS; //fixed
		tmp.X0Ch()=img.X0()*nReS; //fixed
		return tmp;
	}
	else return img;
}



template<typename T>
voxelImageT<T>  resampleMode(const voxelImageT<T>& img, double nReSampleNotSafe)//  TODO to be tested
{
	voxelImageT<T> tmp;
	if (nReSampleNotSafe < .999)  {
		int nReS=1./nReSampleNotSafe+0.5;
		tmp.reset(nReS*img.size3());
		forAllkji_(tmp)
			tmp(i,j,k)=img((0.5+i)/nReS,(0.5+j)/nReS,(0.5+k)/nReS);

		tmp.dxCh()=img.dx()/nReS; //fixed
		tmp.X0Ch()=img.X0()/nReS; //fixed
		return tmp;
	}
	else if (nReSampleNotSafe > 1.001)  {
		int nReS=nReSampleNotSafe+0.5; /// Warning unsigned doesn't work  wTf
		tmp.reset(img.size3()/nReS);
		forAllkji_(tmp)  {
			std::map<T,short> neis;///.  ID-counter
			const T pID = tmp(i,j,k);
			forAllNei(0,nReS-1)
				if(T vj = _nei(img,i*nReS,j*nReS,k*nReS); vj != pID )
					++(neis.insert({vj,0}).first->second);// else ++nSames;
			tmp(i,j,k)=std::max_element(neis.begin(), neis.end(), mapComparer<T>())->first;
		}
		tmp.dxCh()=img.dx()*nReS; //fixed
		tmp.X0Ch()=img.X0()*nReS; //fixed
	return tmp;
	}
	else return img;
}


template<typename T>
voxelImageT<T> resampleMax(const voxelImageT<T>& img, double nReSampleNotSafe)//  TODO to be tested
{
	voxelImageT<T> tmp;
	if (nReSampleNotSafe < .999)  {
		int nReS=1./nReSampleNotSafe+0.5;
		tmp.reset(nReS*img.size3());
		forAllkji_(img)
			tmp(i,j,k)=img((0.5+i)/nReS,(0.5+j)/nReS,(0.5+k)/nReS);

		tmp.dxCh()/=nReS; //fixed
	}
	else if (nReSampleNotSafe > 1.001)  {
		int nReS=nReSampleNotSafe+0.5; /// Warning unsigned doesn't work  wTf
		tmp.reset(img.size3()/nReS);
		forAllkji_(img)  {
			T neiSum=std::numeric_limits<T>::min();
			forAllNei(0,nReS-1)
			{
				neiSum=std::max(neiSum, _nei(img,i*nReS,j*nReS,k*nReS));
			}
			tmp(i,j,k)=neiSum;//(0.5+double(neiSum)/(nReS*nReS*nReS));
		}

		tmp.dxCh()*=nReS;
	}
	return tmp;
}


template<typename T>
void voxelImageT<T>::rotate(char direction)  {// wrong X0
	int n1,n2,n3;

	(std::cout<<" x<->"<<direction<<" ").flush();
	voxelField<T>::getSize(n1,n2,n3);
	if (direction=='z' || direction=='Z')  {
		{
			double X0Tmp=X0_.x;  X0_.x=X0_.z;  X0_.z=X0Tmp;
			double dxTmp=dx_.x;   dx_.x=dx_.z;  dx_.z=dxTmp;
		}
		voxelImageT<T> tmp=*this;
		this->reset(n3,n2,n1,T(0));
		size_t nij =this->nij_;
		OMPFor()
		for (int k=0; k<n3; ++k)
			for (int j=0; j<n2; ++j)
			{
				//for (int i=0; i<n1; ++i) (*this)(k,j,i)=tmp(i,j,k);
				for (int i=1; i<n1 ; i+=2)
				{
					const T& vv0 = tmp(i,j,k);
					const T  vv1 = *(&vv0-1);
					*(&((*this)(k,j,i)=vv0)-nij)=vv1;
				}
				if(n1%2) (*this)(k,j,n1-1)=tmp(n1-1,j,k);
			}
	}
	else if (direction=='y' || direction=='Y')  {
		{
			double X0Tmp=X0_.x;  X0_.x=X0_.y;  X0_.y=X0Tmp;
			double dxTmp=dx_.x;  dx_.x=dx_.y;  dx_.y=dxTmp;
		}

		voxelImageT<T> tmp=*this;
		this->reset(n2,n1,n3,T(0));
		for (int k=0; k<n3; ++k)
			for (int j=0; j<n2; ++j) {
				//for (int i=1; i<n1; ++i) (*this)(j,i,k)=tmp(i,j,k);
				for (int i=1; i<n1 ; i+=2) {
					const T& vv0 = tmp(i,j,k);
					const T  vv1 = *(&vv0-1);
					*(&((*this)(j,i,k)=vv0)-(*this).nnn_.x)=vv1;
				}
				if(n1%2) (*this)(j,n1-1,k)=tmp(n1-1,j,k);
			}

	}
	else if (direction=='-')  {
		std::cout<<" -> flipping image,  x origin will be invalid "<<std::endl;
		voxelImageT<T> tmp=*this;
		for (int k=0; k<n3; ++k)
			for (int j=0; j<n2; ++j)
				for (int i=0; i<n1; ++i)
					(*this)(n1-1-i, j, k)=tmp(i,j,k);
	}
	else
		std::cout<<"\n\nSwapping "<<direction<<" and x directions(!?!), skipping  \n"<<std::endl;

}


template<typename T>
void voxelImageT<T>::PointMedian032(int nAdj0,int nAdj1, T lbl0, T lbl1)  {
	unsigned long nChanges(0);
	int3 nnn = voxelField<T>::size3();
	for (int i=0; i<3; ++i) nnn[i]=nnn[i]+2;

	voxelImageT<T> vxls(*this);
	//vxls.growBox(1);


	forAllkji_1_sum_(vxls, reduction(+:nChanges))  {
		T& vr = (*this)(i,j,k); //(*this)(i-1,j-1,k-1);
		const T vv = (*this)(i,j,k); //(*this)(i-1,j-1,k-1);
		if(lbl0==vv || vv ==lbl1)  {
			int neiSum0=0, neiSum1=0;
			forAllNei(-1,1) {
				T vj=_nei(vxls,i,j,k);
				neiSum0 += vj==lbl0;
				neiSum1 += vj==lbl1;
			}
			T* vp = &vxls(i,j,k);
			T
			vj = vxls.v_i(-1,vp);  neiSum0 += vj==lbl0;  neiSum1 += vj==lbl1;
			vj = vxls.v_i( 1,vp);  neiSum0 += vj==lbl0;  neiSum1 += vj==lbl1;
			vj = vxls.v_j(-1,vp);  neiSum0 += vj==lbl0;  neiSum1 += vj==lbl1;
			vj = vxls.v_j( 1,vp);  neiSum0 += vj==lbl0;  neiSum1 += vj==lbl1;
			vj = vxls.v_k(-1,vp);  neiSum0 += vj==lbl0;  neiSum1 += vj==lbl1;
			vj = vxls.v_k( 1,vp);  neiSum0 += vj==lbl0;  neiSum1 += vj==lbl1;

			//neiSum-=vxls(i,j,k);
			if (neiSum0 >= nAdj0 && neiSum0>neiSum1 && vv==lbl1) {  vr=lbl0; ++nChanges;  }
			else
			if (neiSum1 >= nAdj1 && neiSum1>neiSum0 && vv==lbl0) {  vr=lbl1;  ++nChanges;  }
		}
	}

	std::cout<<"PointMedian032  changed: "<<nChanges<<std::endl;

}



template<typename T>
void FaceMedGrowToFrom(voxelImageT<T>& vImg, T lblTo, T lblFrm, int ndif=0)  {
	unsigned long nChanges(0);
	int3 nnn = vImg.size3();
	for (int i=0; i<3; ++i) nnn[i]=nnn[i]+2;

	voxelImageT<T> vxls=vImg;
	//vxls.growBox(1);

	forAllkji_1_sum_(vxls, reduction(+:nChanges))  {
		const T vv = vImg(i,j,k);////(*this)(i-1,j-1,k-1);
		if(vv ==lblFrm)  {
			int neiSumTo=0, neiSum1=0;

			//short nSames(0);
			std::map<T,short> neis;///.  ID-counter

			T* vp = &vxls(i,j,k);
			T
			vj = vxls.v_i(-1,vp);  neiSumTo += vj==lblTo;  neiSum1 += vj==lblFrm;
			vj = vxls.v_i( 1,vp);  neiSumTo += vj==lblTo;  neiSum1 += vj==lblFrm;
			vj = vxls.v_j(-1,vp);  neiSumTo += vj==lblTo;  neiSum1 += vj==lblFrm;
			vj = vxls.v_j( 1,vp);  neiSumTo += vj==lblTo;  neiSum1 += vj==lblFrm;
			vj = vxls.v_k(-1,vp);  neiSumTo += vj==lblTo;  neiSum1 += vj==lblFrm;
			vj = vxls.v_k( 1,vp);  neiSumTo += vj==lblTo;  neiSum1 += vj==lblFrm;

			if ( neiSumTo && neiSumTo > neiSum1+ndif) {  vImg(i,j,k)=lblTo; ++nChanges;  }
			//else
			//if (neiSum1 >= neiSumTo  && vv==lbl0) {  *vp=lbl1;  ++nChanges;  }
		}
	}

	std::cout<<"FaceMedGrowToFrom  nChanges: "<<nChanges<<std::endl;
}


template<typename T>
size_t  voxelImageT<T>::FaceMedian06(int nAdj0,int nAdj1)  {
	unsigned long long nChanges(0);

	voxelImageT<T> vxls=*this;
	vxls.growBox(1);

	forAllkji_1_sum_(vxls, reduction(+:nChanges))  {
		int neiSum = sum6Nei(vxls, &vxls(i,j,k));

		if (neiSum <= nAdj0 && (*this)(i-1,j-1,k-1))  {
			(*this)(i-1,j-1,k-1)=0;
			++nChanges;
		}
		else if (neiSum >= nAdj1 && !((*this)(i-1,j-1,k-1)))  {
			(*this)(i-1,j-1,k-1)=1;
			++nChanges;
		}
	}

	std::cout<<"FaceMedian06  nChanges: "<<nChanges<<std::endl;
	return nChanges;
}



template<typename T>
voxelImageT<T> medianx(const voxelImageT<T>& vImage)  {
	//unsigned long nChanges(0);
	(std::cout<<"  median ").flush();
	voxelImageT<T> vxls=vImage;
	forAllkji_1_(vImage)  {  const T* vp=&vImage(i,j,k);
		std::array<T,3> vvs={{ *vp,
								vImage.v_i(-1,vp), vImage.v_i( 1,vp)
								}};

		std::nth_element(vvs.begin(),vvs.begin()+1,vvs.end());
		//nChanges+=vxls(i,j,k) != vvs[3];
		vxls(i,j,k) = vvs[1];
	}

	//(std::cout<<nChanges<<", ").flush();
	return vxls;
}


template<typename T>
void voxelImageT<T>::shrinkPore()  {
	voxelImageT<T> vxls=*this;

	forAllkji_1_(vxls)  {
		T* vp = &vxls(i,j,k);
		if (*vp==0 && ( sum6Nei(vxls, vp) ))   (*this)(i,j,k)=1;
	}


	OMPFor()
	for (int k=1; k<vxls.nz()-1; ++k)  {
		for (int j=1; j<vxls.ny()-1; ++j)  {
			//for (int i=0; i<vxls.nx()-1; ++i)//
			{	int i=0;
				if (vxls(i,j,k)==0 && ( vxls(i,j+1,k) || vxls(i,j-1,k) || vxls(i,j,k+1) || vxls(i,j,k-1) ))
					(*this)(i,j,k)=1;

				i=vxls.nx()-1;
				if (vxls(i,j,k)==0 && ( vxls(i,j+1,k) || vxls(i,j-1,k) || vxls(i,j,k+1) || vxls(i,j,k-1) ))
					(*this)(i,j,k)=1;
		   }
		}
	}


	OMPFor()
	for (int k=1; k<vxls.nz()-1; ++k)
		//for (int j=0; j<vxls.ny()-1; ++j)
			for (int i=1; i<vxls.nx()-1; ++i)
			{
				int j=0;
				if (vxls(i,j,k)==0 && ( vxls(i-1,j,k) || vxls(i+1,j,k) || vxls(i,j,k-1) || vxls(i,j,k+1) ))
					(*this)(i,j,k)=1;

				j=vxls.ny()-1;
				if (vxls(i,j,k)==0 && ( vxls(i-1,j,k) || vxls(i+1,j,k) || vxls(i,j,k-1) || vxls(i,j,k+1) ))
					(*this)(i,j,k)=1;
		   }

	//for (int k=0; k<vxls.nz()-1; ++k)
		OMPFor()
		for (int j=1; j<vxls.ny()-1; ++j)
			for (int i=1; i<vxls.nx()-1; ++i)
			{
				int k=0;
				if (vxls(i,j,k)==0 && ( vxls(i-1,j,k) || vxls(i+1,j,k) || vxls(i,j-1,k) || vxls(i,j+1,k) ))
					(*this)(i,j,k)=1;

				k=vxls.nz()-1;
				if (vxls(i,j,k)==0 && ( vxls(i-1,j,k) || vxls(i+1,j,k) || vxls(i,j-1,k) || vxls(i,j+1,k) ))
					(*this)(i,j,k)=1;
		   }

}


template<typename T>
long long mode(voxelImageT<T>& vImg, short minDif, bool verbose)  {
	short nSameSkip=3-(minDif/2);
	long long nChanges = 0;
	voxelImageT<T> vxls=vImg;
	forAllkji_1_sum_(vxls, reduction(+:nChanges))  {
		T* vp = &vxls(i,j,k);
		const T pID = *vp;

		short nSames(0);
		std::map<T,short> neis;///.  ID-counter

		T vj;
		if ( pID != (vj = vxls.v_i(-1,vp)) ) 	 ++(neis.insert({vj,0}).first->second); else ++nSames;
		if ( pID != (vj = vxls.v_i( 1,vp)) ) 	 ++(neis.insert({vj,0}).first->second); else ++nSames;
		if ( pID != (vj = vxls.v_j(-1,vp)) ) 	 ++(neis.insert({vj,0}).first->second); else ++nSames;
		if ( pID != (vj = vxls.v_j( 1,vp)) ) 	 ++(neis.insert({vj,0}).first->second); else ++nSames;
		if ( pID != (vj = vxls.v_k(-1,vp)) ) 	 ++(neis.insert({vj,0}).first->second); else ++nSames;
		if ( pID != (vj = vxls.v_k( 1,vp)) ) 	 ++(neis.insert({vj,0}).first->second); else ++nSames;

		if(nSames<nSameSkip)  {
			auto neitr = std::max_element(neis.begin(), neis.end(), mapComparer<T>());
			if (neitr->second >= nSames+minDif) {
				++nChanges;
				vImg(i,j,k) = neitr->first;
			}
		}
	}
	if(verbose) std::cout<<"  modeMinDif"<<minDif<<"_nChanges:"<< nChanges<<"; ";
	return nChanges;
}

template<typename T>
void voxelImageT<T>::mode(short minDif, bool verbose)  { // ICL version compatibility
	::mode(*this, minDif, verbose);
}

template<typename T>
long long modeNSames(voxelImageT<T>& vImage, const short nSameNei, bool verbose=false)  {
	long long nChanges = 0;
	voxelImageT<T> vxls=vImage;
	forAllkji_1_sum_(vxls, reduction(+:nChanges))  {
		T* vp = &(vxls(i,j,k));
		const T pID = *vp;

		short nSames(0);
		std::map<T,short> neis;///.  ID-counter

		if (T vj = vxls.v_i(-1,vp); vj != pID)  ++(neis.insert({vj,0}).first->second); else ++nSames;
		if (T vj = vxls.v_i( 1,vp); vj != pID)  ++(neis.insert({vj,0}).first->second); else ++nSames;
		if (T vj = vxls.v_j(-1,vp); vj != pID)  ++(neis.insert({vj,0}).first->second); else ++nSames;
		if (T vj = vxls.v_j( 1,vp); vj != pID)  ++(neis.insert({vj,0}).first->second); else ++nSames;
		if (T vj = vxls.v_k(-1,vp); vj != pID)  ++(neis.insert({vj,0}).first->second); else ++nSames;
		if (T vj = vxls.v_k( 1,vp); vj != pID)  ++(neis.insert({vj,0}).first->second); else ++nSames;

		if(nSames<=nSameNei)  {
			 auto neitr = std::max_element(neis.begin(), neis.end(), mapComparer<T>());
			 if (neitr->second > nSames)  {
				++nChanges;
				vImage(i,j,k) = neitr->first;
			 }
		 }
	}
	if(verbose) std::cout<<"  modeNSames("<<nSameNei<<")_nChanges:"<< nChanges<<"; ";
	return nChanges;
}


template<typename T>
void mode26(voxelImageT<T>& vImg, short minDif, bool verbose)  {
	const short minDifp2=minDif/2;
	voxelImageT<T> voxls=vImg;
	long long nChanges = 0;
	forAllkji_1_sum_(voxls, reduction(+:nChanges))  {
		T* vp = &(voxls(i,j,k));
		const T pID = *vp;

		short nSames(0);
		T* neiPID;


		#define neis_countE_m_(j_m,k_m) \
			neiPID = &(voxls(i-1,j_m,k_m));	\
			if (*(  neiPID) == pID  ) ++nSames; \
			if (*(++neiPID) == pID  ) ++nSames; \
			if (*(++neiPID) == pID  ) ++nSames;

		{ neis_countE_m_(j-1,k-1) }
		{ neis_countE_m_(j  ,k-1) }
		{ neis_countE_m_(j+1,k-1) }
		{ neis_countE_m_(j-1,k  ) }
		{ neis_countE_m_(j  ,k  ) }
		{ neis_countE_m_(j+1,k  ) }
		{ neis_countE_m_(j-1,k+1) }
		{ neis_countE_m_(j  ,k+1) }
		{ neis_countE_m_(j+1,k+1) }


		if(nSames+minDifp2<14) {
			std::map<T,short> neis;///.  ID-counter

			#define neis_insert_neiPID_m_(j_m,k_m) \
				neiPID = &(voxls(i-1,j_m,k_m));	\
				if (*(  neiPID) != pID )  ++(neis.insert({*neiPID,0}).first->second); \
				if (*(++neiPID) != pID )  ++(neis.insert({*neiPID,0}).first->second); \
				if (*(++neiPID) != pID )  ++(neis.insert({*neiPID,0}).first->second);

			{ neis_insert_neiPID_m_(j-1,k-1) }
			{ neis_insert_neiPID_m_(j  ,k-1) }
			{ neis_insert_neiPID_m_(j+1,k-1) }
			{ neis_insert_neiPID_m_(j-1,k  ) }
			{ neis_insert_neiPID_m_(j  ,k  ) }
			{ neis_insert_neiPID_m_(j+1,k  ) }
			{ neis_insert_neiPID_m_(j-1,k+1) }
			{ neis_insert_neiPID_m_(j  ,k+1) }
			{ neis_insert_neiPID_m_(j+1,k+1) }

			auto neitr = std::max_element(neis.begin(), neis.end(), mapComparer<T>());
			if (neitr->second > nSames+minDif) {
				vImg(i,j,k) = neitr->first; ++nChanges;
			}
		 }
	}
	if(verbose) std::cout<<"  mode26_"<<minDif<<"_nChanged:"<< nChanges<<"; ";
}



template<typename T>
void voxelImageT<T>::growPore()  {// optimized function, should be further optimized as it is frequently used

	voxelImageT<T> vxls=*this;


	forAllkji_1_(vxls)  {
		if (vxls(i,j,k) && ( !vxls(i-1,j,k) || !vxls(i+1,j,k) ||
							  !vxls(i,j-1,k) || !vxls(i,j,k+1) || !vxls(i,j+1,k) || !vxls(i,j,k-1) ))
			(*this)(i,j,k)=0;
	}


	OMPFor()
	for (int k=1; k<vxls.nz()-1; ++k)
		for (int j=1; j<vxls.ny()-1; ++j) {
			// for (int i=0; i<vxls.nx()-1; ++i)//
				int i=0;
				if (vxls(i,j,k) && ( !vxls(i+1,j,k) || !vxls(i,j-1,k)  || !vxls(i,j,k+1) || !vxls(i,j+1,k) || !vxls(i,j,k-1) ))
					(*this)(i,j,k)=0;

				i=vxls.nx()-1;
				if (vxls(i,j,k) && ( !vxls(i-1,j,k) || !vxls(i,j-1,k)  || !vxls(i,j,k+1) || !vxls(i,j+1,k) || !vxls(i,j,k-1) ))
					(*this)(i,j,k)=0;
		}


	OMPFor()
	for (int k=1; k<vxls.nz()-1; ++k)
		// for (int j=0; j<vxls.ny()-1; ++j)
			for (int i=1; i<vxls.nx()-1; ++i) {
				int j=0;

				if (vxls(i,j,k) && ( !vxls(i-1,j,k) || !vxls(i+1,j,k)  || !vxls(i,j,k+1) || !vxls(i,j+1,k) || !vxls(i,j,k-1) ))
					(*this)(i,j,k)=0;

				j=vxls.ny()-1;
				if (vxls(i,j,k) && ( !vxls(i-1,j,k) || !vxls(i+1,j,k) || !vxls(i,j-1,k)  || !vxls(i,j,k+1) || !vxls(i,j,k-1) ))
					(*this)(i,j,k)=0;
		 }

	// for (int k=0; k<vxls.nz()-1; ++k)
		OMPFor()
		for (int j=1; j<vxls.ny()-1; ++j)
			for (int i=1; i<vxls.nx()-1; ++i) {
				int k=0;

				if (vxls(i,j,k) && ( !vxls(i-1,j,k) || !vxls(i+1,j,k) || !vxls(i,j,k+1) || !vxls(i,j+1,k) || !vxls(i,j-1,k) ))
					(*this)(i,j,k)=0;

				k=vxls.nz()-1;
				if (k>0 && vxls(i,j,k) && ( !vxls(i-1,j,k) || !vxls(i+1,j,k) || !vxls(i,j-1,k) || !vxls(i,j+1,k) || !vxls(i,j,k-1) ))
					(*this)(i,j,k)=0;
			}

}








template<typename T>  void voxelImageT<T>::NOT(const voxelImageT& data2)  {
	forAlliii_((*this))	(*this)(iii)= (*this)(iii) && !data2(iii);
}
template<typename T>  void voxelImageT<T>::AND(const voxelImageT& data2)  {
	forAlliii_((*this))
				(*this)(iii)= (*this)(iii) && data2(iii);
}
template<typename T>  void voxelImageT<T>::OR(const voxelImageT& data2)  {
	forAlliii_((*this))
				(*this)(iii)= (*this)(iii) || data2(iii);
}

template<typename T>  void voxelImageT<T>::XOR(const voxelImageT& data2)  {
	forAlliii_((*this))
		(*this)(iii)= (*this)(iii) != data2(iii);
}

template<typename T>  void voxelImageT<T>::maxEq(const voxelImageT& data2)  {
	forAlliii_((*this))
				(*this)(iii)= max((*this)(iii), data2(iii));
}

template<typename T>  void voxelImageT<T>::minEq(const voxelImageT& data2)  {
	forAlliii_((*this))
				(*this)(iii)= min((*this)(iii), data2(iii));
}


template<typename T>
void voxelImageT<T>::threshold101(T Bgn,T  End)  {
	forAllvp_((*this))  {  T vv = *vp;   *vp= vv<Bgn || End<vv;  }
}

template<typename T,  enable_if_t<std::is_arithmetic<T>::value, int> = 0>
void rescale(voxelImageT<T>& img, T newMin,T  newMax)  {
	T vmin = std::numeric_limits<T>::max();
	T vmax = std::numeric_limits<T>::min();
	int deltaT = newMax - newMin;

	OMPFor(reduction(min:vmin) reduction(max:vmax))
	forAllcp(img) {vmin = std::min(*cp,vmin); vmax = std::max(*cp,vmin); }
	std::cout<<"   vmin:"<<int(vmin)<<"   vmax:"<<int(vmax)<<"  ";
	vmax = std::max(T(vmin+1),vmax);
	forAllvp_(img)  *vp = newMin + (deltaT*(*vp - vmin))/(vmax-vmin);
}




template<typename T>
void voxelImageT<T>::fillHoles(int maxHoleRadius)  { // TODO optimize
	std::cout<<"  filling small isolated parts: "<<std::flush;
	voxelImageT<T> dataTmp=*this;
	std::cout<<"-"<<std::flush;

	dataTmp.shrinkPore(); std::cout<<".";std::cout.flush();
	for (int i=0; i<6; ++i)  {
		dataTmp.growPore();  dataTmp.OR(*this); std::cout<<".";std::cout.flush(); }
	*this=dataTmp;
	std::cout<<"-"<<std::flush;

	dataTmp.growPore();
	for (int i=0; i<4; ++i)  {
		dataTmp.shrinkPore(); dataTmp.AND(*this); std::cout<<".";std::cout.flush();}
	*this=dataTmp;
	std::cout<<"-"<<std::flush;

	if ( maxHoleRadius > 1)  {
		for (int i=0; i<maxHoleRadius; ++i)
			{ dataTmp.shrinkPore(); std::cout<<".";std::cout.flush();}
		for (int i=0; i<maxHoleRadius*4; ++i)
			{ dataTmp.growPore();  dataTmp.OR(*this); std::cout<<".";std::cout.flush();}
		*this=dataTmp;
		std::cout<<"-"<<std::flush;

		for (int i=0; i<maxHoleRadius; ++i)
			{ dataTmp.growPore(); std::cout<<".";std::cout.flush();}
		for (int i=0; i<maxHoleRadius*3; ++i)
			{ dataTmp.shrinkPore(); dataTmp.AND(*this); std::cout<<".";std::cout.flush();}
		*this=dataTmp;
		std::cout<<"-"<<std::flush;
	}
	std::cout<<"."<<std::endl;

}

template<typename T>
void voxelImageT<T>::writeAConnectedPoreVoxel(std::string fnam) const  {
	std::cout<<" finding a connected pore voxel:";
	for (int nShrink=4; nShrink>=0  ; nShrink--)  {
		voxelImageT<T> dataTmp=*this;
		for (int iShrink=1; iShrink<=nShrink ; iShrink++)  {
			dataTmp.shrinkPore();
		}

		std::cout<<" nShrink: "<<nShrink<<std::endl;;
		dataTmp.printInfo();

		for (int k=dataTmp.nz()*3/4-2; k>dataTmp.nz()*1/8+1  ; k--)
		 for (int j=dataTmp.ny()*3/4-2; j>dataTmp.ny()*1/8+1  ; j--)
			for (int i=dataTmp.nx()*3/4-2; i>dataTmp.nx()*1/8+1  ; i--)
			{
				if (dataTmp(i,j,k)==0 &&
					 !(   dataTmp(i+1,j+1,k+1) || dataTmp(i-1,j+1,k+1) || dataTmp(i-1,j-1,k+1) || dataTmp(i+1,j-1,k+1)
						|| dataTmp(i-1,j-1,k-1) || dataTmp(i+1,j-1,k-1) || dataTmp(i-1,j+1,k-1) || dataTmp(i+1,j+1,k-1)
						|| dataTmp(i-1,j,k) || dataTmp(i+1,j,k)
						|| dataTmp(i,j-1,k) || dataTmp(i,j+1,k)
						|| dataTmp(i,j,k-1) || dataTmp(i,j,k+1)
					 )
					)
				{
					std::ofstream of(fnam);
					assert(of);
					of<<(i+0.5)*dx_.x+X0_.x<<" "<<(j+0.5)*dx_.y+X0_.y<<" "<<(k+0.5)*dx_.z+X0_.z<<std::endl;
					of.close();
					std::cout<<" found  ("<<i<<" "<<j<<" "<<k<<") -> "<<double(dataTmp(i,j,k))<<std::endl;
					std::cout<<" found  xyz:  "<<i*dx_.x+X0_.x<<" "<<j*dx_.y+X0_.y<<" "<<k*dx_.z+X0_.z<<std::endl;
					return ;
				}
			}
	}
	std::cout<<" \n			  ----  ERRORR   -----	  \n"<<std::endl;
	std::cout<<" \n\n ----  didn't find  a connected pore voxel   -----  \n\n"<<std::endl;
}

template<typename T>
void replaceRange(voxelImageT<T>& vImage, T minvi, T  maxvi, T midvi)  {
	(std::cout<<"  "<<Tint(minvi)<<":"<<Tint(maxvi)<<"->"<<Tint(midvi)<<"  ").flush();
	const T minv=minvi,   maxv=maxvi,  midv=midvi;
	forAllvp_(vImage)  if (minv<=*vp && *vp<=maxv)  *vp=midv;
}




template<typename T, enable_if_t<!std::is_arithmetic<T>::value, int> = 0 >
void printInfo(const voxelImageT<T>& vImage){}

template<typename T, enable_if_t<std::is_arithmetic<T>::value, int> = 0 >
void printInfo(const voxelImageT<T>& vImage)  {
	int3 nnn=vImage.size3();
	if constexpr (std::is_integral<T>::value)  { // sync "total_porosity" with tests
		unsigned long long nPores=0;
		unsigned long long nTotal=0;
		std::cout<<"\n  Image info, (void==0, valid!="<<Tint(maxT(T))<<"):"<<std::endl;
		OMPFor(reduction(+:nPores) reduction(+:nTotal))
		forAllcp(vImage)  { nPores += (*cp==0); nTotal += (*cp!=maxT(T)); }
		// sync "totalPorosity" with tests
		std::cout << "   totalPorosity: "<< double(nPores)/(double(nnn.x)*nnn.y*nnn.z)  <<"  = "<< nPores<<"/ ("<<nnn.x<<"*"<<nnn.y<<"*"<<nnn.z<<")" << std::endl;
		std::cout << "   validPorosity: "<< double(nPores)/double(nTotal) <<"  = " << nPores<<"/"<<nTotal << std::endl;
		std::cout << "   dx: " << vImage.dx()<<",  X0: "<< vImage.X0() << std::endl;

		int minv=1000000000, maxv=-1000000000; long long avgv=0;
		forAllcp(vImage)  { int val = *cp; minv=std::min(minv,val); maxv=std::max(maxv,val); avgv+=val; }
		std::cout << "   min: "<<minv  << " max: "<<maxv  << " avg: "<<avgv/(double(nnn.x)*nnn.y*nnn.z)<< std::endl;
	}
	else {
		double minv=1e64, maxv=-1e64, avgv=0.;
		OMPFor(reduction(min:minv) reduction(max:maxv) reduction(+:avgv))
		forAllcp(vImage)  { double val = toScalar(*cp); minv=std::min(minv,val); maxv=std::max(maxv,val); avgv+=val;

			ensure(std::isfinite(avgv),"overflow or ???",-1);
		}
		std::cout << "   min: "<<minv  << " max: "<<maxv  << " avg: "<<avgv/(double(nnn.x)*nnn.y*nnn.z)<< std::endl;
	}

}

template<typename T>  void voxelImageT<T>::printInfo() const   {  ::printInfo(*this);  }


template<typename T>
double voxelImageT<T>::volFraction(T vv1, T vv2) const {
	unsigned long long nPores=0;
	int3 nnn=(*this).size3();
	OMPFor(reduction(+:nPores))
	forAllcp((*this))	nPores += ( vv1<=*cp && *cp<=vv2 );

	return double(nPores)/(double(nnn.x)*nnn.y*nnn.z);
}


template<typename T>
voxelImageT<T> median(const voxelImageT<T>& vImage)  {
	//unsigned long nChanges(0);
	(std::cout<<"  median ").flush();
	voxelImageT<T> vxls=vImage;
	forAllkji_1_(vImage)  {  const T* vp=&vImage(i,j,k);
		std::array<T,7> vvs={{ *vp,
								vImage.v_i(-1,vp), vImage.v_i( 1,vp),
								vImage.v_j(-1,vp), vImage.v_j( 1,vp),
								vImage.v_k(-1,vp), vImage.v_k( 1,vp)
								}};

		std::nth_element(vvs.begin(),vvs.begin()+3,vvs.end());
		//nChanges+=vxls(i,j,k) != vvs[3];
		vxls(i,j,k) = vvs[3];
	}

	//(std::cout<<nChanges<<", ").flush();
	return vxls;
}


template<typename T>
void circleOut(voxelImageT<T>& vImage, int X0,int Y0,int R, char dir='z', T outVal=maxT(T))  {
	int rlim = R*R;
	if (dir=='z')  {
	 forAllkji_(vImage)
		  if((i-X0)*(i-X0)+(j-Y0)*(j-Y0)>rlim)
		    vImage(i,j,k)=outVal;
	}
	else if (dir=='x')  {
	 for (int k=0; k<vImage.nz(); k++)
	  for (int j=0; j<vImage.ny(); ++j)
		if((j-X0)*(j-X0)+(k-Y0)*(k-Y0)>rlim)
	      for (int i=0; i<vImage.nx(); ++i)
		    vImage(i,j,k)=outVal;
	} else alert("Bad circleOut direction "+_s(dir),-1);
}



template<typename T>
void maskWriteFraction(voxelImageT<T>& vImage, std::string maskname, std::string fnam, unsigned char maskvv, T minIelm, T maxIelm) {
	//  TODO to be tested
	voxelImageT<unsigned char> mask(maskname);
	T maxvv = std::min(maxIelm, accumulate(vImage,(std::max<T>)));//(const T& (*)(const T&, const T&))(std::max<T>)
	std::cout<<"  maxvv:"<<maxvv<<std::endl;
	std::vector<int> nMasked(maxvv+3,0);
	std::vector<int> nNotmsk(maxvv+3,0);

	forAllkji_(vImage)  {	T vv=vImage(i,j,k);
		if(minIelm<=vv && vv<=maxIelm)  {
			if(mask(i,j,k)==maskvv)		++nMasked[vv];
			else                       ++nNotmsk[vv];
		}
	}
	std::cout<<" Mask Info:"<<std::endl;
	mask.printInfo();
	//mask.write("dumpMask.mhd");
	std::ofstream of(fnam); ensure(of);
	if(of)  {
		std::cout<<"  Writting "<<fnam<<std::endl;
		for(T i=minIelm; i<=maxvv; ++i) of<<double(nMasked[i])/(nMasked[i]+nNotmsk[i]+1e-38)<<std::endl;
	}
}


template<typename T>
void mapToFrom(voxelImageT<T>& vImage, const voxelImageT<T>& vimg2, T vmin, T vmax, double scale=0, double shift=0.5)  { // no interpolation
	int3 N2=vImage.size3();
	size_t count=0;
	int3 N1;
	int3 Dn((vImage.X0()-vimg2.X0())/vImage.dx()+0.5);
	N1.x=std::max(-Dn.x,0);
	N1.y=std::max(-Dn.y,0);
	N1.z=std::max(-Dn.z,0);
	N2.x=std::min(vimg2.nx()-Dn.x, N2.x);
	N2.y=std::min(vimg2.ny()-Dn.y, N2.y);
	N2.z=std::min(vimg2.nz()-Dn.z, N2.z);

	std::cout << " mapping bounds: "<<Dn<<" to "<<N2<<std::endl;
	OMPFor(reduction(+:count))
	for (int k=N1.z; k<N2.z; ++k)
	for (int j=N1.y; j<N2.y; ++j)
	for (int i=N1.x; i<N2.x; ++i)  {
		T&  vv = vImage(i, j, k);
		if(vmin<=vv && vv<=vmax)  {
			vv = vimg2(i+Dn.x,j+Dn.y,k+Dn.z);
			++count;
		}
		if(scale>1e-16)  vv = vv*scale+shift;
	}
	std::cout << "  N Changed: "<<count<<",  "<<100.*count/(double(N2.x-N1.x)*(N2.y-N1.y)*(N2.z-N1.z))<<"%"<<std::endl;
}

template<typename T>
void mapToFrom(voxelImageT<T>& vImage, const voxelImageT<T>& vimg2)  { // no interpolation
	int3 N2=vImage.size3();
	int3 N1;
	int3 Dn((vImage.X0()-vimg2.X0())/vImage.dx()+0.5);
	N1.x=std::max(-Dn.x,0);
	N1.y=std::max(-Dn.y,0);
	N1.z=std::max(-Dn.z,0);
	N2.x=std::min(vimg2.nx()-Dn.x, N2.x);
	N2.y=std::min(vimg2.ny()-Dn.y, N2.y);
	N2.z=std::min(vimg2.nz()-Dn.z, N2.z);

	std::cout << " mapping bounds: "<<Dn<<" to "<<N2<<std::endl;
	OMPFor()
	for (int k=N1.z; k<N2.z; ++k)
	for (int j=N1.y; j<N2.y; ++j)
	for (int i=N1.x; i<N2.x; ++i)  {
			vImage(i, j, k) = vimg2(i+Dn.x,j+Dn.y,k+Dn.z);
	}
}





template<typename T>
  typename std::enable_if<std::is_integral<T>::value,bool>::type
  operat( voxelImageT<T>& img1, char opr, const voxelImageT<T>& img2, std::stringstream& ins)  {

	double shift(0.);  if(!std::isalpha(opr)) ins>>shift;
	int shiftI=shift+0.5;

	ensure(img1.nx()==img2.nx()," image sizes do not match "+_s(img1.size3())+" != "+_s(img2.size3()),-1);
	ensure(img1.size3()==img2.size3()," image sizes do not match "+_s(img1.size3())+" != "+_s(img2.size3()));

	if(shift) std::cout<<" "<<opr<<":shift:"<<shift<<std::endl; else std::cout<<" ";

	switch (opr) {
		case '+':  forAlliii_(img1) img1(iii)=std::max(Tint(minT(T)),std::min(Tint(img1(iii))+(img2(iii)-shiftI),Tint(maxT(T))));	break;
		case '-':  forAlliii_(img1) img1(iii)=std::max(Tint(minT(T)),std::min(Tint(img1(iii))-(img2(iii)-shiftI),Tint(maxT(T))));	break; //! Analyse needs this
		case '&':  forAlliii_(img1) { img1(iii)= img1(iii)&img2(iii); }  break;
		case '|':  forAlliii_(img1) { img1(iii)= img1(iii)|img2(iii); }  break;
		case '%':  forAlliii_(img1) { img1(iii)= img1(iii)%img2(iii); }  break;
		case '*':  shift+=1.; forAlliii_(img1) { img1(iii)=std::min(Tint(double(img1(iii))*img2(iii)*shift),Tint(maxT(T))); }  break;
		case '/':  shift+=1.; forAlliii_(img1) { img1(iii)=std::min(Tint(double(img1(iii))/img2(iii)*shift),Tint(maxT(T))); }  break;
		case 'm':  // map masked
		{
			int3 X0(0,0,0);
			int       msk1Bgni=0,      msk1Endi=maxT(T),     msk2Bgni=0,      msk2Endi=maxT(T);
			ins>>X0>> msk1Bgni       >>msk1Endi        >>msk2Bgni       >>msk2Endi;
			T msk1Bgn=msk1Bgni,msk1End=msk1Endi, msk2Bgn=msk2Bgni,msk2End=msk2Endi;
																	ensure(msk1Bgn<=msk1End);  ensure(X0.x+img2.nx()<=img1.nx());
			std::cout<<" X0:"<<X0<<" msk1:"<<int(msk1Bgn)<<".."<<int(msk1End)<<"   msk2:"<<msk2Bgni<<".."<<msk2Endi<<std::endl;
			forAllkji_(img2) {
				const T v2=img2(i,j,k);
				if(msk2Bgn<=v2 && v2<=msk2End) {
					T& v1=img1(X0.x+i,X0.y+j,X0.z+k);
					if(msk1Bgn<=v1 && v1<=msk1End)  v1=v2;	}	}
			break;
		}
		default:   std::cout<<" !!!operation "<<opr<<" voxelImage not supported!!! ";
			return 1;
	}
	return 0;
}

template<typename T>
  typename std::enable_if<std::is_integral<T>::value,bool>::type
  operat( voxelImageT<T>& vImg, char opr, std::string img2Nam, std::stringstream& ins)  {

	if(img2Nam.empty())  {
	 std::cout<<"  image ="<<opr<<"image ";
	 switch (opr) {
		case '!':
			forAllvp_(vImg) { *vp= !(*vp); }  break;
		case '~':
			forAllvp_(vImg) { *vp= ~(*vp); }  break;
		case '-':
			forAllvp_(vImg) { *vp= -(*vp); }  break;
		default:   std::cout<<" !!!operation "<<opr<<" not supported!!! ";
	 }
	}else{
	 if (!isalpha(img2Nam[0])&&!isalpha(img2Nam.back()))
	 {
		Tint   ii; ii=strTo<Tint>(img2Nam);
		double dd; dd=strTo<double>(img2Nam);
		std::cout<<"  image "<<opr<<"= "<< dd<<" ";
		switch (opr) {
			case '=':
				forAllvp_(vImg) { (*vp)= ii; }	break;
			case '+':
				forAllvp_(vImg) { (*vp)+= ii; }	break;
			case '-':
				forAllvp_(vImg) { (*vp)-= ii; }	break;
			case '&':
				forAllvp_(vImg) { (*vp)= (*vp)&ii; }	break;
			case '|':
				forAllvp_(vImg) { (*vp)= (*vp)|ii; }	break;
			case '%':
				forAllvp_(vImg) { (*vp)= (*vp)%ii; }	break;
			case '*':
				if(dd>1) forAllvp_(vImg) { (*vp)= std::min(Tint((*vp)*dd),Tint(maxT(T))); }
				else    { forAllvp_(vImg) { (*vp)*= dd; } }	break;
			case '/':
				forAllvp_(vImg) { (*vp)/= dd; }	break;
			case 'b':
				forAllvp_(vImg) { (*vp)= std::max(ii,Tint(*vp)); }	break;
			case 'e':
				forAllvp_(vImg) { (*vp)= std::min(ii,Tint(*vp)); }	break;
			default:   std::cout<<" !!!operation "<<opr<<" "<<img2Nam<<" not supported!!! ";
		}
	 }else{
		std::cout<<"  image  "<<opr<<"= "<<img2Nam<<" ";
		_dbgetOrReadImg(T, img2, img2Nam);
		operat(vImg,opr,img2,ins);
	 }
	}


	return true;
}



template<typename T> typename std::enable_if<!std::is_integral<T>::value,int>::type
  operat( voxelImageT<T>& vImg, char opr, std::string img2Nam, std::stringstream & ins)  {

	double shift(0.);  ins>>shift;
	if(img2Nam.empty())  {
	 std::cout<<"  image ="<<opr<<"image ";
	 switch (opr) {
		case '-':
			forAllvp_(vImg) { *vp *= -1; }  break;
		default:   std::cout<<"  Image  "<<opr<<"= !!!not supported!!! ";
	 }
	}
	else  {
	 if (!isalpha(img2Nam[0])&&!isalpha(img2Nam.back()))
	 {
		Tint   ii; ii=strTo<Tint>(img2Nam);
		double dd; dd=strTo<double>(img2Nam);
		std::cout<<"  image "<<opr<<"= "<< ii<<" ";
		switch (opr) {
			case '=':
				forAllvp_(vImg) { (*vp) =ii; }  break;
			case '+':
				forAllvp_(vImg) { (*vp)+=ii; }  break;
			case '-':
				forAllvp_(vImg) { (*vp)-=ii; }  break;
			case '*':
				forAllvp_(vImg) { (*vp)*=dd; }  break;
			case '/':
				forAllvp_(vImg) { (*vp)/=dd; }  break;
			//case 'b':
				//forAllvp_(vImg) (*vp)=std::max(ii,(*vp));	break;
			//case 'e':
				//forAllvp_(vImg) (*vp)=std::min(ii,(*vp));	break;
			default:   std::cout<<" !!!operation "<<opr<<" not supported!!! ";
		}
	 }else{
		std::cout<<"  image  "<<opr<<"= "<<img2Nam<<" ";
		_dbgetOrReadImg(T, img2,img2Nam);

		switch (opr) {
			case '+':
				forAlliii_(vImg) { vImg(iii)=vImg(iii)+img2(iii); }	break;
			case '-':
				//forAlliii_(vImg) vImg(iii)=min(max(int(vImg(iii))+128-img2(iii),0),maxT(T));	break;
				forAlliii_(vImg) { vImg(iii)=vImg(iii)-img2(iii); }	break; //! Analyse needs this
			case '*':
				forAlliii_(vImg) { vImg(iii)=vImg(iii)*img2(iii); }	break;
			//case '/':
				//forAlliii_(vImg) vImg(iii)/=img2(iii);	break;
			default:   std::cout<<" !!!operation "<<opr<<" not supported!!! ";
		}
	 }
	}

	return 0;
}
