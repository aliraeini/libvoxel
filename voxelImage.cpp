/*-------------------------------------------------------------------------*\

This file is part of libvoxel, a C++ template library for handelling 3D images.

Developed by:
 - Ali Q Raeini (2010-2022)

You can redistribute this code and/or modify this code under the
terms of the GNU General Public License (GPL) as published by the
Free Software Foundation, either version 3 of the License, or (at
your option) any later version. see <http://www.gnu.org/licenses/>.

\*-------------------------------------------------------------------------*/

#include <memory>
#include <sstream>
#include "globals.h"  // ensure...
#include "InputFile.h"

#include "voxelImageI.h"
#include "voxelEndian.h"
#include "voxelRegions.h"

// backward compatibility with ICL-2021 versions TODO: delete
#if MS_VERSION < 2
#include "vxlPro.cpp"
#endif

using namespace std; //cin cout endl string stringstream  istream istringstream regex*


int maxNz = 50000;


namespace MCTProcessing{
	template<typename T> VxlFuncs<T> namedProcesses();
}

template<typename T>
class  voxelplugins {
 public:
	typedef bool(*ProcessP)( stringstream&  inputs, voxelImageT<T>& vImg);

	VxlFuncs<T> key_funs;
	const VxlFuncs<T>& operator()() const { return key_funs; }
	voxelplugins() {
		using namespace MCTProcessing;
		key_funs = MCTProcessing::namedProcesses<T>();
	};


	int vxProcess(const InputFile& inp, voxelImageT<T>& img, string nam="") const  { // nam is ignored here
		if(inp.data().size()>2) std::cout<<std::endl;
		for(const auto& ky:inp.data())  {
			auto paer = key_funs.find(ky.first);
			if (paer!=key_funs.end())  {
				(cout<<" "<<ky.first<<": ").flush();
				stringstream ss(ky.second);
				(*(paer->second))(ss, img);
				if(inp.data().size()>2) cout<<endl;
			}
			else  {
				if(ky.first!="end") { cout<<"  stopped executing "+inp.fileName()+" before \""+ky.first+"\"  :/ ";
												return -1; }
				break;
			}
		}
		return 0;
	};
	int vxProcess(const string&  keystr, voxelImageT<T>& img, string nam)	const {  return  vxProcess(InputFile(keystr,nam,false), img);   }
	/*int process(      istream& keyins, voxelImageT<T>& img, string nam)	{  return  process(InputFile(keyins,nam,false), img);
		//while (true)  {
			//std::streampos begLine = keyins.tellg();
			//string ky;  keyins>>ky;
			//if (keyins.fail()) 	{cout<<"  @"<<keyins.tellg()<<"  "<<nam<<" done."<<endl;  break; }
			//else if (ky[0]=='{' || ky[0]=='}') { keyins.seekg(int(keyins.tellg())-ky.size()-1); continue; }
			//else if (ky[0]=='#' || ky[0]=='\'' || ky[0]=='/' || ky[0]=='%')  keyins.ignore(10000,'\n');
			//else  {
				//auto paer = key_funs.find(ky);
				//if (paer!=key_funs.end())  {
					//(cout<<" "<<ky<<": ").flush();
					//stringstream ss;  if(keyins.peek()!='\n') keyins.get (*(ss.rdbuf()));
					//(*(paer->second))(ss,img);
					//cout<<endl;
				//}
				//else  {
					//cout<<"  stopped processing "<<nam<<" before \""<<ky<<"\" :/ "<<endl;
					//keyins.clear(); keyins.seekg(begLine);
					//return -1;
				//}
		//}	}
		//return 0;
	};*/

};

template<class InpT, typename T>  //! run voxel plugins
 int vxlProcess(const InpT& ins, voxelImageT<T>& img, string nam) {
	return voxelplugins<T>().vxProcess(ins,img,nam);  }

template<class InpT, typename First, typename... Rest>
 int vxlProcess(const InpT& ins, voxelImageTBase* imgPtr, string nam)  { //! detect type and run voxel plugins
	if(auto img = dynamic_cast<voxelImageT<First>*>(imgPtr))
		return vxlProcess<InpT,First>(ins,*img,nam);
	else if(sizeof...(Rest))
		return vxlProcess<InpT,Rest...>(ins, imgPtr, nam);
	cout<<"Unknown image type."<<endl;
	return -1;
}

//TODO: break down and parallelize compilation
template int vxlProcess<InputFile, SupportedVoxTyps>(const InputFile& inp, voxelImageTBase* imgPtr, string nam);
template int vxlProcess<string, SupportedVoxTyps>(string const& ins, voxelImageTBase* imgPtr, string nam);


string VxlKeysHelp(string keyname, string subkey)  {
	//! Query and print MCTProcessing keyword usage messages

	VxlFuncs<unsigned char> key_funs = voxelplugins<unsigned char>()();

	voxelImage vImg;
	stringstream keys;
	if(keyname.size())  {
		auto paer = key_funs.find(keyname);
		if (paer!=key_funs.end())  {
			stringstream ss(subkey.empty()?  "?" : "? "+subkey);
			try                         {  (*(paer->second))(ss, vImg); }
			catch (std::exception &exc) {  std::cerr <<keyname<<" KeyHelp not implemented:" << exc.what() << endl; }
			catch (...)                 {  std::cerr <<keyname<<" KeyHelp not implemented:" << endl; }
			return ss.str();
		}
		else
			cout<<" Error: no such keyword "<<keyname<<endl;
		keys<<"//!-*- C -*- keywords:\n";
		for(const auto& proc:key_funs) 	keys<<proc.first<<"\n";
		keys<<" Error: no such keyword "<<keyname<<"\n\n";
	}
	else  {
		std::vector<std::pair<string,VxlFunc<unsigned char>>> keyfuns(key_funs.begin(), key_funs.end());
		std::sort(keyfuns.begin(), keyfuns.end());
		for(const auto& proc:keyfuns)  if(proc.first.size()>1)  {
			stringstream ss("?");
			try                         {  (*(proc.second))(ss, vImg); }
			catch (std::exception &exc) {  std::cerr <<proc.first<<" KeyHelp not implemented:" << exc.what() << endl; }
			catch (...)                 {  std::cerr <<proc.first<<" KeyHelp not implemented:" << endl; }
			keys<<"  "<< std::left<<std::setw(15)<<proc.first+":"<<" "<<replaceFromTo(ss.str(),"\n","\n                  ")<<"\n";
		}
	}
	return keys.str();
}


template<typename T>
void voxelImageT<T>::readFromHeader(const string& hdrNam, int procesKeys)  {
	//! read image from file header, format detected based on image extension

	if (hdrNam.empty() || hdrNam=="NO_READ")  return;

	std::ifstream fil{hdrNam};  ensure(fil,"Cannot open header file, "+hdrNam,-1);
	//! read image from file header, format detected based on image extension
	auto& vImg=*this; string fnam;

	int3 nnn(0,0,0);
	string BinaryData="XXX", flipSigByt="False";
	bool X0read=false, dxread=false, autoUnit=true; //auto unit only applies to .mhd format
	double unit_=1.;
	int nSkipBytes(0);
	if (hasExt(hdrNam,".mhd")) {
		(cout<<" "<<hdrNam<<": ").flush();
		while (true)  {
			std::streampos begLine = fil.tellg();
			string ky, tmp;   fil>>ky>>tmp;
			stringstream ss;  if(fil.peek()!='\n') fil.get (*(ss.rdbuf()));
			if (fil.fail()) break;
			if (ky=="ObjectType")  {  ss>> tmp;  if (tmp != "Image") cout<<"  Warning: ObjectType != Image :="<<tmp<<endl;	}
			else if (ky=="NDims")  {  ss>> tmp;  if (tmp != "3"    ) cout<<"  Warning: NDims != 3 :="<<tmp<<endl;	}
			else if (ky=="ElementType")  { ss>> tmp;  if ((tmp != "MET_UCHAR") && (sizeof(T)==1)) cout<<"  Warning: ElementType != MET_UCHAR :="<<tmp<<endl; 	}
			else if (ky=="Offset")       { ss>> vImg.X0_;   cout<<"  X0: "<<  vImg.X0_<<",  ";  X0read=true; }
			else if (ky=="ElementSize"
			      || ky=="ElementSpacing")  { ss>> vImg.dx_;  cout<<"  dX: "<<vImg.dx_<<",  ";  dxread=true;	}
			else if (ky=="DimSize")                              {  ss>> nnn;  nnn.z=std::min(nnn.z,maxNz);  cout<<"  Nxyz: "<<nnn<<",  ";	}
			else if (ky=="ElementDataFile")  {  if (fnam.empty()) ss>> fnam;
			                                  if (size_t is=hdrNam.find_last_of("\\/"); is<hdrNam.size() && fnam[0]!='/' && fnam[1]!=':')
			                                      fnam=hdrNam.substr(0,is+1)+fnam;
			                                   cout<<"  Img: "<<fnam<<",	"; }
			else if (ky=="BinaryData")  {  ss>> BinaryData;     cout<<"  BinaryData: "<<BinaryData<<"	"<<endl; }
			else if (ky=="Unit")        {  ss>> unit_;  autoUnit=false;   cout<<"  Unit, OneMeter: "<<unit_<<endl; 	}
			else if (ky=="HeaderSize")  {  ss>> nSkipBytes;         cout<<"  Ski pHeaderSize: "<<nSkipBytes<<endl;	}
			else if (ky=="OutputFormat" || ky=="DefaultImageFormat" )  {  if(tmp=="=") ss>> tmp;  cout<<"  OutputFormat: "<<tmp<<", suffix:"<<imgExt(tmp)<<"	"<<endl; }///. sets suffix+format
			else if (ky=="BinaryDataByteOrderMSB" || ky=="ElementByteOrderMSB")  {  ss>> flipSigByt; }
			else if (ky!="CompressedData" &&  ky!="CompressedDataSize" &&  ky!="TransformMatrix" &&
					 ky!="ElementNumberOfChannels" && ky!="CenterOfRotation" && ky!="AnatomicalOrientation" && ky!="AnatomicalOrientation")  {
				fil.clear();  fil.seekg(begLine);
				(cout<<"; ").flush();
				break;
			}
		}
		cout<<endl;
	}
	#ifdef TIFLIB
	else if (hasExt(hdrNam,".tif"))  {  readTif(vImg, hdrNam);  return;  }
	#endif
	else if (hasExt(hdrNam,".am"))	{
		fnam=hdrNam;
		procesKeys=0;
	}
	else if (hasExt(hdrNam,".raw.gz") || hasExt(hdrNam,".raw") || hasExt(hdrNam,".dat"))  { // detect size and voxel size from image name.
		string
		data=replaceFromTo(replaceFromTo(replaceFromTo(replaceFromTo(replaceFromTo(
									hdrNam,".gz$",""), ".raw$",""), ".dat$",""),"__","\n"),"_"," ");
		data=replaceFromTo(replaceFromTo(replaceFromTo(data,"voxel",""),"size",""),"um","\n");
		data=regex_replace(data,regex("( [0-9][0-9]*)c"), " $1 $1 $1 ", regex_constants::format_first_only);
		data=regex_replace(data,regex("( [0-9][0-9]*)[ x]*([0-9][0-9]*)[ x]*([0-9][0-9]* )"),
		                                        "\n   reset Nd0 $1 $2 $3 ", regex_constants::format_first_only);
		data=regex_replace(data,regex("^[^\n]*\n"), "", regex_constants::format_first_only);
		data=regex_replace(data,regex("\n|($)"),"\n   read "+hdrNam+"\n", regex_constants::format_first_only);
		for(auto&da:data)  { if(da=='p') da='.'; else if(da=='\n') break; }
		cout<<"  Keywords: {\n"<<data<<"  }"<<endl;
		vxlProcess(data,vImg,hdrNam);
		procesKeys=0;
	}
	else if (hasExt(hdrNam,"_header"))  {
		cout<<" (depricated) _header:"<<hdrNam<<","<<endl;

		char tmpc;
		for (int i=0; i<8; ++i)   fil>>tmpc, cout<<tmpc;  //ignore the first 8 characters (ascii 3uc)

		if (hasExt(hdrNam,"_header"))  fnam=hdrNam.substr(0,hdrNam.size()-7);
		fil>>nnn >> vImg.dx_ >>	vImg.X0_ ;
		cout<<"\n Nxyz: "<<nnn<<"    dX: "<< vImg.dx_<<"   X0: "<< vImg.X0_ <<" um"<< endl;
		ensure(fil,"incomplete/bad header name", -1);
	}
	else  alert("Unknown (header) file type: "+hdrNam,-1); // exit

	if(nnn.z) vImg.reset(nnn);
	int readingImage=0;
	if( !fnam.empty() && fnam!="NO_READ" && procesKeys!=2)  {
	  if (hasExt(fnam,".tif")) {
			dbl3 dx=vImg.dx_, X0=vImg.X0_;
			readingImage = vImg.readBin(fnam);
			if(X0read) vImg.X0_=X0;
			if(dxread) vImg.dx_=dx;
	  }
	  else if ((hasExt(fnam,".raw") && BinaryData!="False") || BinaryData=="True")   {
			readingImage = vImg.readBin(fnam, nSkipBytes);
	  }
	  else if (hasExt(fnam,".am"))    {
			int RLECompressed;
			dbl3 dx=vImg.dx_, X0=vImg.X0_;
			getAmiraHeaderSize(fnam, nnn,vImg.dx_,vImg.X0_,nSkipBytes,RLECompressed);
			readingImage = vImg.readBin(fnam, nSkipBytes);
			if(X0read) vImg.X0_=X0;
			if(dxread) vImg.dx_=dx;
	  }
	  else if (hasExt(fnam,".raw.gz")) {
			readingImage = vImg.readBin(fnam);
	  }
	  else   {
		std::ifstream in(fnam);  assert(in);
		if(nSkipBytes) in.ignore(nSkipBytes);
		vImg.voxelField<T>::readAscii(in);
	  }
	}
	ensure(readingImage==0, "cannot read image "+fnam,-1);

	if(flipSigByt=="True") {
		cout<<"  flipEndian "<<endl;
		flipEndian(vImg);	}

	if(autoUnit  && vImg.dx_[0]>0.02)	{ //&& dxread
		cout<<"   dx="<<vImg.dx_[0]<<"(>0.02 -> assuming unit is um), ";
		unit_ = 1e-6;
	}
	vImg.dx_*=unit_;
	vImg.X0_*=unit_;
	if(abs(unit_-1.)>epsT(float)) cout<<"  unit= "<<unit_<<" => dx= "<<vImg.dx_<<", X0= "<<vImg.X0_<<endl;


	if (procesKeys) voxelplugins<T>().vxProcess(InputFile(fil,hdrNam),vImg);
	cout<<"."<<endl;
}


std::unique_ptr<voxelImageTBase> readImage(string hdrNam,	int procesKeys)  {
	//! read or create image
	using namespace std;
	(cout<<"voxelImage \""<<hdrNam<<"\": ").flush();
	if (hasExt(hdrNam,".am"))  {
		string vtype = getAmiraDataType(hdrNam);
		cout<<"reading '"<<vtype<<"'s from .am file"<<endl;

		if (vtype=="int")       return make_unique<voxelImageT<int>>(hdrNam,0);
#ifdef _ExtraVxlTypes
		if (vtype=="short")     return make_unique<voxelImageT<short>>(hdrNam,0);
#endif
		if (vtype=="ushort")    return make_unique<voxelImageT<unsigned short>>(hdrNam,0);
		if (vtype=="byte")      return make_unique<voxelImageT<unsigned char>>(hdrNam,0);

		alert("data type "+vtype+" not supported, when reading "+hdrNam, -1);
	}

	#ifdef TIFLIB
	if (hasExt(hdrNam,".tif"))  return readTif(hdrNam);
	#endif

	string typ;
	std::ifstream fil(hdrNam); // header file
	if(!fil)
	{
		ensure(hdrNam.size()<4 || hdrNam[hdrNam.size()-4]!='.', "can not open header file '"+hdrNam+"', pwd: "+getpwd(), -1);
		typ = hdrNam; hdrNam="NO_READ";
	}
	else if (hasExt(hdrNam,".mhd")) {
		while (true)  {
			string ky;  fil>>ky;
			stringstream ss;
			if(fil.peek()!='\n') fil.get (*(ss.rdbuf()));
			if (fil.fail()) {  cout<<"\n\n\nWarning: readImage, 'ElementType =' not set in "<<hdrNam<<endl; break; }
			if (ky == "ElementType")  {  ss >> typ >> typ;  break; }
		}
	}
	fil.close();

	if (typ=="MET_UCHAR")        return make_unique<voxelImageT<unsigned char>>(hdrNam, procesKeys);
	if (typ=="MET_USHORT")       return make_unique<voxelImageT<unsigned short>>(hdrNam, procesKeys);
	if (typ=="MET_INT")          return make_unique<voxelImageT<int>>           (hdrNam, procesKeys);
	if (typ=="MET_FLOAT")        return make_unique<voxelImageT<float>>         (hdrNam, procesKeys);
#ifdef _ExtraVxlTypes
	if (typ=="MET_CHAR")         return make_unique<voxelImageT<char>>          (hdrNam, procesKeys);
	if (typ=="MET_SHORT")        return make_unique<voxelImageT<short>>         (hdrNam, procesKeys);
	if (typ=="MET_UINT")         return make_unique<voxelImageT<unsigned int>>  (hdrNam, procesKeys);
	if (typ=="MET_DOUBLE")       return make_unique<voxelImageT<double>>        (hdrNam, procesKeys);
	//if (typ=="MET_FLOAT_ARRAY")  return make_unique<voxelImageT<float3>>        (hdrNam, procesKeys);
	//if (typ=="MET_DOUBLE_ARRAY") return make_unique<voxelImageT<dbl3>>          (hdrNam, procesKeys);
#endif //_ExtraVxlTypes
	return                              make_unique<voxelImage>(hdrNam, procesKeys);

}
