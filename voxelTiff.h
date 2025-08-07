#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <typeinfo>
#include <typeindex>
#include <unordered_map>
#include <string.h>
#include <iostream>

#include "tiffio.h"
#include "voxelImage.h"

//https://stackoverflow.com/questions/24059421/adding-custom-tags-to-a-tiff-file, didn't work properly, also a too dirty solution

template<typename T>
int tifDataType(T) {
	const std::unordered_map<std::type_index, int> tifTypes
	{
		{std::type_index(typeid(unsigned char)), SAMPLEFORMAT_UINT},
		{std::type_index(typeid(char)),          SAMPLEFORMAT_INT},
		{std::type_index(typeid(int)),           SAMPLEFORMAT_INT},
		{std::type_index(typeid(unsigned int)),  SAMPLEFORMAT_UINT},
		{std::type_index(typeid(short)),         SAMPLEFORMAT_INT},
		{std::type_index(typeid(unsigned short)), SAMPLEFORMAT_UINT},
		{std::type_index(typeid(float)),         SAMPLEFORMAT_IEEEFP},
		{std::type_index(typeid(double)),        SAMPLEFORMAT_IEEEFP},
	};
	return tifTypes.at(std::type_index(typeid(T)));
}


inline void getTifTags(dbl3& X0_, dbl3& dx_, TIFF *tif) {
	float dx=0., dy=0.;
	TIFFGetField(tif, TIFFTAG_XRESOLUTION, &dx);
	TIFFGetField(tif, TIFFTAG_YRESOLUTION, &dy);
	dx_={dx,dy,dy};

	float X0=0., Y0=0., Z0=-2.1e30;
	TIFFGetField(tif, TIFFTAG_XPOSITION, &X0);
	TIFFGetField(tif, TIFFTAG_YPOSITION, &Y0);

	uint32_t kFrst=0;
	TIFFGetField(tif, TIFFTAG_IMAGEDEPTH, &kFrst);
	Z0=double(kFrst)*dx;

	X0_={X0,Y0,Z0};

	char * info;
	if(TIFFGetField(tif, TIFFTAG_IMAGEDESCRIPTION, &info)) {
		std::istringstream instrim(info);
		while(instrim.good())  {
			std::string str;
			instrim>>str;
			if(str=="dx")	instrim >> dx_;
			else if(str=="X0")	instrim >> X0_;
		}
	}
	if(dx_.x<1e-16) (std::cout<<"\n !!! Warning: dx read from tif seems invalid  !!!\n").flush();
}

inline void setTifTags(const dbl3& X0_, const dbl3& dx_, TIFF *tif) {
	// negative numbers not supported in TIFFTAG_XPOSITION :-( for X0, using TIFFTAG_IMAGEDESCRIPTION instead
	float dx=dx_[0], dy=dx_[1];
	TIFFSetField(tif, TIFFTAG_XRESOLUTION, dx);
	TIFFSetField(tif, TIFFTAG_YRESOLUTION, dy);
	//TIFFSetField(tif, TIFFTAG_ZRESOLUTION, dz);, dz=dx_[2]

	std::ostringstream ostrim;
	ostrim<<" dx "<<dx_<<" ";
	if(mag(X0_)>1e-16) ostrim<<" X0 "<<X0_<<" ";
	(std::cout<<" tag: \""<<ostrim.str()<<"\" ").flush();
	TIFFSetField(tif, TIFFTAG_IMAGEDESCRIPTION, ostrim.str().c_str());

	TIFFSetField(tif, TIFFTAG_SOFTWARE, "voxelImage+libtiff");
}


template<typename T>
int readTif( voxelField<T>&  aa, std::string fnam )
{
	(std::cout<<  " reading tif file "<<fnam<<" ").flush();

	uint32_t nx, ny;

	TIFF *tif = (TIFF *) NULL;
	tif = TIFFOpen(fnam.c_str(), "r");	 if (tif == NULL)	{ alert("Cannot open tif"); return -1; }

	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &nx);
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &ny);
	int npages=TIFFNumberOfDirectories(tif);
	if (npages>1.1*maxNz+2)  npages = maxNz;

	(std::cout<<"size: "<<aa.size3()<<" * "<<sizeof(T)).flush();
	if (voxelImageT<T>* vxls = dynamic_cast<voxelImageT<T>*>(&aa)) {
		getTifTags(vxls->X0Ch(),vxls->dxCh(),tif);
		(std::cout<<",  X0:"<<vxls->X0()<<",  dx:"<<vxls->dx()).flush(); }

	aa.reset(nx,ny,npages,T(0.));


	for (int pn=0;pn<npages;++pn) {
		TIFFReadEncodedStrip( tif, static_cast<tstrip_t>(0), static_cast<void *>(&aa(0,0,pn)), static_cast<tsize_t>(nx * ny) * sizeof(T) );
		//for (row = 0; row < ny; row++) {		if (TIFFReadScanline(tif, &aa(0,row,pn), row, 0) < 0)		break;	}
		TIFFReadDirectory(tif);
	}
  	TIFFClose(tif);

	std::cout<<  " ."<<std::endl;
	return 0;
}


template<typename T>
int writeTif(const voxelField<T>&  aa, std::string fnam) {

	uint32_t nx=aa.nx(), ny=aa.ny();
	int pn=0, npages=aa.nz();
	int smplfrmt = tifDataType(T());

	TIFF * tif = (TIFF *) NULL;
	tif = TIFFOpen(fnam.c_str(), "w8");		if (tif == NULL)	return -2;

	const voxelImageT<T>* vxls = dynamic_cast<const voxelImageT<T>*>(&aa);
	if(vxls) setTifTags(vxls->X0(),vxls->dx(),tif);
	else std::cout<<"dxXo not set"<<std::endl;

	for(pn=0; pn<npages; ++pn) {

		TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, nx);
		TIFFSetField(tif, TIFFTAG_IMAGELENGTH, ny);

		TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8*sizeof(T));
		TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
		TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, smplfrmt);
		TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK); // this is needed by paraview!

		TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW); // For best compression, choose COMPRESSION_PACKBITS and manually compress as 7z. Other opts:COMPRESSION_NONE  COMPRESSION_DEFLATE COMPRESSION_PACKBITS

		TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, ny);

		TIFFSetField(tif, TIFFTAG_PAGENUMBER, pn, npages);

		T aOrig=aa(nx/2,ny/2,pn);  // Warn if image modified in libtiff
		TIFFWriteEncodedStrip(tif, 0, const_cast<T *>(&aa(0,0,pn)), nx*ny*sizeof(T));
		if(aOrig!=aa(nx/2,ny/2,pn)) std::cout<<"Warning image modified in libtiff"<<std::endl;

		TIFFWriteDirectory(tif);
	}

	TIFFClose(tif);
	return (0);
}


template<typename T>
int writeTif(const voxelField<T>&  aa, std::string fnam, int iStart,int iEnd , int jStart,int jEnd , int kStart,int kEnd ) {

	voxelImageT<T>  bb;
	bb.reset({iEnd-iStart, jEnd-jStart, kEnd-kStart});
	const voxelImageT<T>* vxls = dynamic_cast<const voxelImageT<T>*>(&aa);
	if(vxls) bb.setFrom(*vxls, iStart, jStart, kStart);
	else     bb.voxelField<T>::setFrom(aa, iStart, jStart, kStart);

	return writeTif(bb, fnam);

	// Deleted approach complained about LZW compression.
	// Note although the deflate compression algorithm overall is a better choice,
	// for now we have to stick to LZW for portability
	// (specifically because Avizo does not support inflate algorithm).
}


inline std::unique_ptr<voxelImageTBase>  readTif(std::string fnam) {
	TIFF *tif = (TIFF *) NULL;
	tif = TIFFOpen(fnam.c_str(), "r");
	if (tif == NULL)  return std::unique_ptr<voxelImageTBase>();
	uint16_t frmt = SAMPLEFORMAT_UINT,   nbits = 8;
	TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &frmt);
	TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &nbits);
	TIFFClose(tif);

	std::cout<<"frmt: "<<frmt<<" : "<<nbits<<std::endl;
	switch (frmt) {
		case SAMPLEFORMAT_UINT:
		default:
			if (nbits==8)  return std::make_unique<voxelImageT<unsigned char>>(fnam);
	#ifndef _VoxBasic8
			if (nbits==16) return std::make_unique<voxelImageT<unsigned short>>(fnam);
#ifdef _ExtraVxlTypes
			if (nbits==32) return std::make_unique<voxelImageT<unsigned int>>(fnam);
#endif
		  break;
		case SAMPLEFORMAT_INT:
#ifdef _ExtraVxlTypes
			if (nbits==16) return std::make_unique<voxelImageT<short>>(fnam);
			if (nbits==8)  return std::make_unique<voxelImageT<char>>(fnam);
#endif
			if (nbits==32) return std::make_unique<voxelImageT<int>>(fnam);
		  break;
		case SAMPLEFORMAT_IEEEFP:
#ifdef _ExtraVxlTypes
			if (nbits==64) return std::make_unique<voxelImageT<double>>(fnam);
#endif
			if (nbits==32) return std::make_unique<voxelImageT<float>>(fnam);
			break;
	#endif
		}
	return nullptr;
}
