#include <memory>
#include <sstream>
#include "voxelImage.h"
#include "voxelImageI.h"
#include "shapeToVoxel.h"
#include "globals.h"  // ensure...
#include "voxelRegions.h"

								namespace MCTProcessing _begins_


template<typename T> ErC ignore( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("{ section to ignore }");
	return 0;
}

template<typename T> ErC fillHoles( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("maxHoleSize // fill small isolated features");
	unsigned int maxHoleSize;
	ins>>maxHoleSize;

	cout<<"  fillHoles: eliminating isolated rocks/pores; maxHoleSize:" <<maxHoleSize<<" (default is 2) "<<endl;
	vImg.fillHoles(maxHoleSize);

	vImg.FaceMedian06(1,5);
	//vImg.FaceMedian07(2,5);
	//vImg.FaceMedian07(2,5);
	return 0;
}

template<typename T> ErC info( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("// print porosity...");
	vImg.printInfo();
	return 0;
}

template<typename T> ErC selectPore( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("thresholdBegin ThresholdEnd;// segment the image");
	cout<<"  Converting to binary (0 and 1):";
	unsigned int  thresholdMin=0,thresholdMax=0;
	ins>>thresholdMin;
	ins>>thresholdMax;

	(cout<<"  pore (=0) <- ["<<int(thresholdMin)<<" "<<int(thresholdMax)<<"] ").flush();
	vImg.threshold101(thresholdMin,std::min(unsigned(maxT(T)),thresholdMax));
	return 0;
}

template<typename T,  enable_if_t<std::is_arithmetic<T>::value, int> = 0>
ErC rescale( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("range_min range_max;// rescale values to be within range");
	(cout<<"  rescaling voxel values to [ ").flush();
	unsigned int  thresholdMin=0,thresholdMax=0;
	ins>>thresholdMin;
	ins>>thresholdMax;

	(cout<<thresholdMin<<", "<<thresholdMax<<" ]    ").flush();
	rescale(vImg,T(thresholdMin),T(thresholdMax));
	(cout<<".").flush();
	return 0;
}

template<typename T> ErC growPore( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("valToGrow(0/1) Algorithm(f) ...;// grow labels to adjacent voxels");
	cout<<"  growing voxels: ";
	int vvalTogrow(0); ins>>vvalTogrow;
	char alg('f'); ins>>alg;

	ensure(vImg.nnn_.z>1, "nz=1, consider  growLabel instead",-1);
	if (!ins.good()) cout<<"\n\n didn't grow anything, usage growPore 0 f f f\n\n"<<endl;
	while (ins.good()) {  // loop while extraction from file is possible
		if (alg=='f')  {
		cout<<"   "<<vvalTogrow<<" "<<alg<<" ";

			if(vvalTogrow==0)
				vImg.growPore();
			else if (vvalTogrow==1)
				vImg.shrinkPore();
			else
			{
				std::cerr<<"growing is only implemented for binary images: "<<
				"selected voxel value to grow is "<<vvalTogrow << ", which is not acceptable"<<endl;
				return false;
			}
		} else {
			std::cerr<<"selected growing algorithm: "<<alg<<
			" the only implemented algorithm is f which stands for faceGrowing"<<endl;
			return false;
		}

		ins>>vvalTogrow;
		ins>>alg;
	}
	cout<<"."<<endl;
	return 0;
}


template<typename T> ErC resampleMean( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("by_resampleFactor // assigning mean val");
	double nResample=1;
	(ins>>nResample, cout<<__FUNCTION__<<" factor: "<<nResample<<" ").flush();
	vImg = resampleMean(vImg,nResample);
	return 0;
}


template<typename T> ErC resampleMax( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("by_resampleFactor // assigning max val");
	double nResample=1;
	(ins>>nResample, cout<<__FUNCTION__<<" factor: "<<nResample<<" ").flush();
	vImg = resampleMax(vImg,nResample);
	return 0;
}

template<typename T> ErC resliceZ( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("nResample(-1)// reslize z to make dz same as dx and dy");
	double nResample=-1;
	(ins>>nResample, cout<<__FUNCTION__<<" factor: "<<nResample<<" ").flush();
	vImg = resliceZ(vImg,nResample);
	return 0;
}

template<typename T> ErC resampleMode( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("nResample //  assigning mode val");
	double nResample=1;
	(ins>>nResample, cout<<__FUNCTION__<<" factor: "<<nResample<<" ").flush();
	vImg = resampleMode(vImg,nResample);
	return 0;
}

template<typename T> ErC redirect( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("direction(y/z) // flip x with y or z axes");
	char axs;
	ins>>axs;
	(cout<<axs<<", swapping x and "<<axs<<" axes").flush();

	vImg.rotate(axs);
	cout<<endl;
	return 0;
}

template<typename T> ErC replaceRange( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("FromVal ToValue NewValue ...");

	Tint  thresholdMin(0),thresholdMax(0), value(0);
	while (ins >> thresholdMin >> thresholdMax >> value) {
	cout<<"  Replacing range  ["<<thresholdMin<<"  "<<thresholdMax<<"] with "<<value<<", ";
	replaceRange(vImg,T(thresholdMin),T(thresholdMax),T(value));
	}
	(cout<<".").flush();
	return 0;
}


template<typename T> ErC cropD( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("cropBegin(0 0 0) cropEnd(nx ny nz)");
	int3 cropBegin(0,0,0), cropEnd=vImg.size3();
	int nLayers(0); int value(1);
	ins>>cropBegin>>cropEnd >> nLayers >> value;
	(cout<<" "<<cropBegin<<" -- "<<cropEnd<<" ").flush();
	if (nLayers) { cout<<"  + "<<nLayers<<" layers of "<<value<<" "<<endl; }
	vImg.cropD(cropBegin,cropEnd,nLayers,value,true);
	return 0;
}

template<typename T> ErC cropf( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("beginFraction(0 0 0) endFraction(0.5 0.5 0.5)");
	dbl3 bgn(0,0,0), end(1,1,1);
	int nLayers(0); int value(1);
	ins>>bgn>>end >> nLayers >> value;
	(cout<<" "<<bgn<<" -- "<<end<<" ").flush();
	if (nLayers) { cout<<"  + "<<nLayers<<" layers of "<<value<<" "<<endl; }
	vImg.cropD(bgn*dbl3(vImg.size3())+0.5,end*dbl3(vImg.size3())+0.5,nLayers,value,true);
	return 0;
}

template<typename T> ErC write( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("outputImageName.raw/mhd/tif/am/.raw.gz");
	string outName("dump.tif");    ins >> outName;
	vImg.write(outName);
	(cout<<".").flush();
	return 0;
}

template<typename T> ErC write8bit( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("outputImageName_8bit.raw/mhd/tif/am/.raw.gz");
	string outName("dump.tif");     ins >> outName;
	double minv=-0.5, maxv=maxT(T);

	{	T mxl = 0;
		OMPragma("omp parallel for reduction(max:mxl)")
		forAllcp(vImg) if(*cp>=0) mxl = max(mxl,*cp);
		if (double(mxl)<255.5) maxv=255.;
	}

	ins>>minv>>maxv;
	double delv=255.499999999/(maxv-minv);
	(cout<<minv<<" "<<maxv).flush();
	voxelImageT<unsigned char> voxels(vImg.size3(),vImg.dx(),vImg.X0(),255);
	forAlliii_(voxels) voxels(iii)=std::max(0,std::min(255,int(delv*(vImg(iii)-minv))));
	voxels.write(outName);
	(cout<<".").flush();
	return 0;
}

template<typename T> ErC read( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("ImageToRead.mhd/.am/.tif");
	int3 nnn = vImg.size3();
	int processHdr=1;  string fnam;   ins>>fnam>>processHdr;
	cout<<"  reading from  image "<<fnam<<endl;
	if(fnam.size()>4)  {
		if ((nnn[2] && (hasExt(fnam,".raw.gz") || hasExt(fnam,".raw"))) || hasExt(fnam,".tif") )  {
			vImg.reset(nnn,T(0));
			vImg.readBin(fnam);
		}
		else vImg.readFromHeader(fnam,processHdr);
	}
	return 0;
}

template<typename T> ErC readAtZ( stringstream& ins, voxelImageT<T>& vImg)  {// used to stitch images
	KeyHint("ImageToReadAndReplacePreviousFromZ.mhd/.am iSlic");
	int3 nnn = vImg.size3();
	size_t iSlic=0;	string fnam;	ins>>fnam>>iSlic;
	cout<<"  reading from  image "<<fnam<<", assigning to slices after "<<iSlic<<endl;
	voxelImageT<T> img(fnam);
	ensure(img.nx()==nnn.x);	ensure(img.ny()==nnn.y);
	std::copy(img.begin(),img.end(),vImg.begin()+iSlic*nnn[0]*nnn[1]);
	return 0;
}



template<typename T> ErC medianFilter( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("nIterations");
	int nIterations(1);  ins >> nIterations;
	(cout<<"  median Filter, nIterations: "<<nIterations).flush();
	vImg.growBox(2);
	for (int i=0; i<nIterations; ++i)  vImg=median(vImg);
	vImg.shrinkBox(2);
	(cout<<".").flush();
	return 0;
}

template<typename T> ErC modeFilter( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("nIterations(1) nMinDif(2) kernelVol(6_or_26)  ");
	int nIterations=1, nMinDif=1, kernelVol=6;   ins >> nIterations >> nMinDif >> kernelVol;
	(cout<<"  mode Filter, nIterations: "<<nIterations<<"  nMinDif"<<nMinDif<<"  kernelVol:"<<kernelVol).flush();
	vImg.growBox(2);
	for (int i=0; i<nIterations; ++i)
		if (kernelVol<10) mode(vImg,nMinDif,true);
		else              mode26(vImg, nMinDif,true);
	vImg.shrinkBox(2);
	(cout<<".").flush();
	return 0;
}

template<typename T> ErC medianX( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("nIterations// only applied in x direction to reduce compressed file size");
	int nIterations=1;   ins >> nIterations;
	(cout<<"  median Filter, nIterations: "<<nIterations).flush();
	for (int i=0; i<nIterations; ++i)
		vImg=medianx(vImg);
	(cout<<".").flush();
	return 0;
}


template<typename T> ErC FaceMedian06( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("nIterations(1),  nAdj0(2), nAdj1(4)");
	int nAdj0(2), nAdj1(4),  nIterations(1);     ins >>nIterations >>nAdj0 >>nAdj1;
	(cout<<"  FaceMedian06: "<<nIterations<<"   "<<nAdj0<<" "<<nAdj1<<"    ").flush();
	vImg.growBox(2);
	for (int i=0; i<nIterations; ++i) vImg.FaceMedian06(nAdj0,nAdj1);
	vImg.shrinkBox(2);
	(cout<<".").flush();
	return 0;
}



template<typename T> ErC PointMedian032( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("nItrs(1),  nAdj(11), lbl0(0), lbl1(1);\n// nAdj is the threshold count of adjacent voxels to change the voxel");
	int nItrs(1),  nAdj(11), lbl0(0), lbl1(1);
	ins >> nItrs >> nAdj >> lbl0 >> lbl1;
	(cout<<"  PointMedian032, "<<" nItrs:"<<nItrs<< "; nAdjThreshold "<<nAdj<<"  lbl0:"<<lbl0<<"  lbl1;"<<lbl1<<"s \n  PointMedian032 is only applied to the labels  lbl0 and  lbl1").flush();
	//vImg.growBox(2);

	for (int i=0; i<nItrs; ++i)  vImg.PointMedian032(nAdj,nAdj,lbl0,lbl1);

	//vImg.shrinkBox(2);
	(cout<<".").flush();
	return 0;
}


template<typename T> ErC faceMedNgrowToFrom( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("nItrs(2),  lblTo(0), lblFrm(1), ndif(-3)");
	int nItrs(2),  ndif(-3);
	Tint lblTo(0), lblFrm(1);
	ins  >> nItrs >> lblTo >> lblFrm >> ndif;
	(cout<<"{ "<<" nItrs:"<<nItrs<<"; "<<lblFrm<<" --> "<<lblTo<< "; ndif: "<<ndif<<";  ").flush();

	vImg.growBox(2); cout<<endl;
	for (int i=0; i<nItrs; ++i) FaceMedGrowToFrom(vImg,T(lblTo),T(lblFrm),ndif);
	vImg.shrinkBox(2);

	(cout<<"};\n").flush();
	return 0;
}

template<typename T> ErC delense032( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("nItrs(2) lbl0(0) lbl1(1) nAdj0(10) nAdj1(6)");
	int nItrs(2),  nAdj0(10),  nAdj1(6);    Tint lbl0(0), lbl1(1);
	ins >> nItrs >> lbl0 >> lbl1 >> nAdj0 >> nAdj1;
	(cout<<"{ "<<" nItrs:"<<nItrs<<"; lbls: "<<lbl0<<" "<<lbl1<< "; nAdjThresholds: "<<nAdj0<<" "<<nAdj1<<";  ").flush();

	vImg.growBox(2); cout<<endl;
	voxelImageT<T> vimgo=vImg;
	for (int i=0; i<nItrs; ++i)   vImg.PointMedian032(25,nAdj1,lbl0,lbl1);
	FaceMedGrowToFrom(vImg,T(lbl1),T(lbl0),1);
	FaceMedGrowToFrom(vImg,T(lbl0),T(lbl1),-1);
	for (int i=0; i<2*nItrs; ++i) { vImg.PointMedian032(nAdj0,25,lbl0,lbl1);	FaceMedGrowToFrom(vImg,T(lbl0),T(lbl1),-1); }
	FaceMedGrowToFrom(vImg,T(lbl0),T(lbl1),-3);
	FaceMedGrowToFrom(vImg,T(lbl0),T(lbl1),-1);
	FaceMedGrowToFrom(vImg,T(lbl0),T(lbl1),-1);

	FaceMedGrowToFrom(vimgo,T(lbl1),T(lbl0),2);//41 51 -> lbl1
	FaceMedGrowToFrom(vimgo,T(lbl1),T(lbl0),2);//41 51 -> lbl1
	forAlliii_(vimgo) if(vimgo(iii)==lbl1) vImg(iii)=lbl1;
	vImg.shrinkBox(2);

	(cout<<"};\n").flush();
	return 0;
}


template<typename T> ErC circleOut( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("Axis(x/y/z) X0(=N/2) Y0(=N/2) R(=N/2) outVal(Max)");

	char d='z';  ins >> d;
	int i = std::max<int>(d-'x',0);
	int X0(vImg.size3()[(i+1)%3]/2), Y0(vImg.size3()[(i+2)%3]/2);
	int R((X0+Y0)/2);
	Tint outVal=maxT(T);

	ins >> X0 >> Y0 >> R >> outVal;
	(cout<<"  circleOut: dir="<<d<<",  X0="<<X0 <<"  Y0="<<Y0  <<"  R="<<R<<"  out="<<outVal ).flush();

	circleOut(vImg,X0,Y0,R,d,T(outVal));

	(cout<<".").flush();
	return 0;
}


template<typename T> ErC maskWriteFraction( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("not implemented");
	int maskvv(2);
	Tint minIelm(1), maxIelm=std::numeric_limits<T>::max();
	string maskname, outName("maskWriteFraction.txt");
	ins >> maskname >> outName >> maskvv >> minIelm >> maxIelm;
	(cout<<"  maskWriteFraction:  mask:"<<maskname <<"  outName:"<<outName<<"  maskvv:"<<maskvv  <<"  minIelm:"<<minIelm<<"  maxIelm:"<<maxIelm ).flush();

	//maskWriteFraction(vImg,maskname,outName,maskvv,minIelm,maxIelm);

	(cout<<"Error: not implemented.").flush();
	return 0;
}


template<typename T> ErC Offset( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("offset(0.,0.,0.)");
	dbl3 offset(0.);  ins >> offset;
	(cout<<"  Offset:"<<offset<<" " ).flush();
	vImg.X0Ch()=offset;
	(cout<<".").flush();
	return 0;
}


template<typename T>  ErC keepLargest0( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint(" // sets smaller isolated regions (of value 0) to 254, computationally expensive");
	keepLargest0(vImg); //! CtrlF:isolated=254
	(cout<<".").flush();
	return 0;
}

template<typename T>  ErC growLabel( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("vvalue(255)  nIters(0) ");
	int  vv(255), nIters(0);  ins >> vv >> nIters;
	(cout<<"  growLabel: "<<vv<<" x"<<nIters ).flush();

	for (int i=0; i<=nIters; ++i)  vImg.growLabel(vv);

	(cout<<".").flush();
	return 0;
}

template<typename T>  ErC reset( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("param Value // param can be N X0 dX V VN or NdXX0");
	string param;
	ins >>param;
	while(param.size())  {
		if(param=="maxNz")  { 	ins>>maxNz;	 cout<<"maxNz:"<<maxNz<<" "; } // stop reading after maxNz layers
		else if(param=="N")  { 	int3 N{0,0,0}; 	ins>>N; 	vImg.reset(N,T{0}); cout<<"N:"<<N<<" "; 	}
		else if(param[0]=='X')  { 	dbl3 X0; 	ins>>X0; 	vImg.X0Ch()=X0; cout<<"X0:"<<X0<<" "; 	}
		else if(param[0]=='d')  { 	dbl3 dx(1.,-2e9,1.); 	ins>>dx;
			if(dx[1]<-1e9) { dx[1]=dx[0]; dx[2]=dx[1]; }
			vImg.dxCh()=dx; cout<<"dX:"<<dx<<" "; 	}
		else if(param[0]=='s')  { 	double scale; 	ins>>scale;
			vImg.dxCh()*=scale; vImg.X0Ch()*=scale; cout<<" dx*="<<scale<<" "; 	}
		else if(param[0]=='V') { // VNdX
		 	int3 N=vImg.size3();
			Tint vv(0.); 	ins>>vv>>N;
			vImg.reset(N,vv);
			if(ins && param=="VNdX") { param="dX"; continue;}
		}
		else if(param[0]=='N' && param[1]=='d') // NdX or Nd0
		{ 	dbl3 dx(1.,-2e9,1.), X0(0.,0.,0.); int3 N=vImg.size3(); 	ins>>N>>dx>>X0;
			if(dx[1]<-1e9) { dx[1]=dx[0]; dx[2]=dx[0]; }
			vImg.reset(N); cout<<"N:"<<N<<" ";  vImg.dxCh()=dx; cout<<"dX:"<<dx<<" "; 	vImg.X0Ch()=X0; cout<<"X0:"<<X0<<" "; 	}
		else cout<<"reset does not support "<<param<<endl;
		param="";
		ins >>param;
	}
	(cout<<".").flush();
	return 0;
}


template<typename T>  ErC operat( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("operation(+-^...) [img2Nam/number]  [shift]");
	string opr=" ", img2Nam;    ins>>opr>>img2Nam;
	operat(vImg,opr[0],img2Nam,ins);
	return 0;
}

template<typename T>
ErC mapFrom( stringstream& ins, voxelImageT<T>& vImg)  {
	KeyHint("image2name  minv maxv");
	Tint minv(0), maxv(maxT(T)); double scale=0, shift(0.5-T(0.5));
	string image2name;
	ins>>image2name>>minv>>maxv;
	ensure(maxv>=minv);
	cout<<"\n{  mapping from image "<<image2name<<", assigning to values originally in range: ["<<minv<<" "<<maxv<<"], by "; if(scale>1e-16) { cout<<shift<<"+"<<scale<<"*"; } cout<<image2name<<endl;
	voxelImageT<T> image2(image2name);

	mapToFrom(vImg,image2,T(minv),T(maxv), scale, shift);

	cout<<" } //mapFrom "<<endl;
	return 0;
}


template<typename T>
VxlFuncs<T> namedProcesses() requires(sizeof(T)<=2) {
	auto key_funs = VxlFuncs<T>{
		{  ""                  , ignore },
		{  ";"                 , ignore }, // TODO delete
		{  "skip"              , ignore  },
		{  "fillHoles"         , fillHoles  },
		{  "reset"             , reset  },
		{  "info"              , info  },
		{  "rescale"           , rescale  },
		{  "pore"              , selectPore  },
		{  "threshold"         , selectPore  },
		{  "threshold101"      , selectPore  },
		{  "Offset"            , Offset  },
		{  "redirect"          , redirect  },
		{  "direction"         , redirect  },
		{  "crop"              , cropD  },
		{  "cropD"             , cropD  },
		{  "cropf"             , cropf  },
		{  "growPore"          , growPore  },
		{  "resample"          , resampleMean  },
		{  "resampleMean"      , resampleMean  },
		{  "resampleMax"       , resampleMax  },
		{  "resampleMode"      , resampleMode  },
		{  "resliceZ"          , resliceZ  },
		{  "rangeTo"           , replaceRange  },
		{  "replaceRange"      , replaceRange  },
		{  "write"             , write  },
		{  "write8bit"         , write8bit  },
		{  "read"              , read  },
		{  "readAtZ"           , readAtZ  },
		{  "modeFilter"        , modeFilter  },
		{  "medianFilter"      , medianFilter  },
		{  "medianX"           , medianX  },
		{  "FaceMedian06"      , FaceMedian06  },
		{  "PointMedian032"    , PointMedian032  },
		{  "faceMedNgrowToFrom", faceMedNgrowToFrom  },
		{  "delense032"        , delense032  },
		{  "circleOut"         , circleOut  },
		{  "growLabel"         , growLabel  },
		{  "keepLargest0"      , keepLargest0  },
		{  "maskWriteFraction" , maskWriteFraction  },
		{  "mapFrom"           , mapFrom  },
		{  "Paint"             , Paint  },
		{  "PaintAdd"          , PaintAdd  },
		{  "PaintBefore"       , PaintBefore  },
		{  "PaintAfter"        , PaintAfter  },
		{  "PaintAddBefore"    , PaintAddBefore  },
		{  "PaintAddAfter"     , PaintAddAfter  },
		{  "operation"         , operat },
		{  "operat"            , operat }
	};

	return key_funs;
}


template<typename T> VxlFuncs<T> namedProcesses() requires(sizeof(T)>=3) {
	auto key_funs = VxlFuncs<T>{
		{  ""          , ignore},// ProcessP can be removed if using g++
		{  ";"         , ignore }, // TODO delete
		{  "skip"	     , ignore },
		{  "Offset"    , Offset },
		{  "redirect"  , redirect },
		{  "direction" , redirect },
		{  "write"     , write },
		{  "read"      , read },
		{  "circleOut" , circleOut }
		};

	return key_funs;
}


template VxlFuncs<unsigned char>  namedProcesses();
template VxlFuncs<unsigned short> namedProcesses();
template VxlFuncs<int>           namedProcesses();
template VxlFuncs<float>          namedProcesses();

								_end_of_(namespace MCTProcessing)
