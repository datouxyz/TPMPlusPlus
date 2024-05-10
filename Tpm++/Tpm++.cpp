// Tpm++.cpp: 定义应用程序的入口点。
//

#include "Tpm++.h"
#include "TissueProblityMap.h"
int GetNumPhyCpu();
#include "TPMInternalFuncs.h"
using namespace std;
#include "struct2x/json/decoder.h"
#include "struct2x/json/encoder.h"
#include <fstream>
// enum TPMOutMask
// {
// AffineMat=0x00000001,FieldInverse=0x00000002,FieldForward=0x00000004,BiasCorrected=0x00000008,
// BiasField=0x00000010,TissImported=0x00000020,TissNative =0x00000040,TissModulate=0x00000080,
// TissUnModulate =0x00000100
// };
std::vector<std::string> SplitString(const std::string& s, char delimiter)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}
sTissues DefaultTissues()
{
	sTissues tss=
	{sTissue(1, true, false, false, false), sTissue(1, true, false, false, false),
		sTissue(2, true, false, false, false), sTissue(3, true, false, false, false),
		sTissue(4, true, false, false, false), sTissue(2, false, false, false, false) };
	return tss;
}




int main(int argc, char *argv[])
{
	if (argc < 5)
	{
		cout << "Usage: " << argv[0] << " -in inputimage -tpm tpmimage -out outFile -json TpmParam " << endl;		
		 cout<<"Out file will has a postfix 'outimage_FieldInverse.nii, outimage_AffineMat.txt,etc...'"<<endl;
		return 1;
	}
	sChannels Channels(1);
	sTpmInputs  TpmInputs;
	ItkImageVec InputImages;
	sitk::Image TpmImg;
	string outimgpath;
	bool bAffineOnly=false;
	// check params for -in -tpm -out -outaffine in a loop
	for(int iarg = 1;iarg<argc;iarg++)
	{
		if (strcmp(argv[iarg], "-in") == 0)
		{
			if (iarg + 1 >= argc)
			{
				cout << "missing input image" << endl;
				return 1;
			}
			// load images until next '-'param
			while (iarg + 1 < argc && argv[iarg + 1][0] != '-')
			{
				sitk::Image img = sitk::ReadImage(argv[iarg+1]);
				InputImages.push_back(img);
				iarg++;
			}
		}
		else if (strcmp(argv[iarg], "-tpm") == 0)
		{
			if (iarg + 1 >= argc)
			{
				cout << "missing tpm image" << endl;
				return 1;
			}
			TpmImg = sitk::ReadImage(argv[iarg+1],sitk::sitkFloat64);
			iarg++;
		}
		else if (strcmp(argv[iarg], "-out") == 0)
		{
			if (iarg + 1 >= argc)
			{
				cout << "missing out outPath" << endl;
				return 1;
			}
			outimgpath = argv[iarg+1];
			iarg++;
		}
		else if (strcmp(argv[iarg], "-json") == 0)
		{
			if (iarg + 1 >= argc)
			{
				cout << "missing json params "<< endl;
				return 1;
			}		
			string jsonstr = LoadStringFromFile(argv[iarg+1]);			
			bool bDecode = struct2x::JSONDecoder(jsonstr.c_str()) >> TpmInputs;			
			iarg++;
		}
		else if (strcmp(argv[iarg], "-AffineOnly") == 0)
		{
			if (iarg + 1 >= argc)
			{
				cout << "missing AffineOnly value"<< endl;
				bAffineOnly = false;
			}
			else
			bAffineOnly = (argv[iarg+1] == string("1")) ? true:false;
		}
	}
	
	if(InputImages.size()==0)
	{	
		cout << "missing input image" << endl;
		return 1;
	}
	if(TpmImg.GetSize().size()==0)
	{
		cout << "missing tpm image" << endl;
		return 1;
	}
	if(outimgpath.size()==0)
	{
		cout << "missing out outimgpath" << endl;
		return 1;
	}	
	

	Eigen::initParallel();
	int ncpu = GetNumPhyCpu();
	Eigen::setNbThreads(ncpu);
	omp_set_num_threads(ncpu);
	while (InputImages.size() > TpmInputs.m_Channels.size())
	{
		TpmInputs.m_Channels.push_back(sChannel());
	}
	if (TpmImg.GetPixelID() != sitk::sitkFloat64)
	{
		ErrOutput(u8"Warn",u8"TPM Will be cast to Float 64(8byte)");
		TpmImg = sitk::Cast(TpmImg,sitk::sitkFloat64);		
	}

	sTPMPriors tpmPriors;
	tpmPriors.LoadPriors(TpmImg);
	
	for (int ichan = 0; ichan < TpmInputs.m_Channels.size(); ichan++)
	{
		auto &chan = TpmInputs.m_Channels[ichan];
		auto & inimg = InputImages[ichan];
		if (inimg.GetPixelID() != sitk::sitkFloat32)
		{
			ErrOutput( u8"Warn",u8"Input Image will be cast to Float 32(4byte)");
			inimg = sitk::Cast(inimg,sitk::sitkFloat32);
		}
		chan.m_Img = inimg;
	}

	MAT_TYPE Affine = TpmInputs.m_InitAffine;
	
	ErrOutput(u8"Affine:", u8"AffinePass Start");
	
	Affine = AffinePass(TpmInputs.m_Channels[0].m_Img, TpmInputs.m_nsample, (TpmInputs.m_fwhm + 1) * 16, tpmPriors, Affine, TpmInputs.m_ptype);
	Affine = AffinePass(TpmInputs.m_Channels[0].m_Img, TpmInputs.m_nsample, TpmInputs.m_fwhm, tpmPriors, Affine, TpmInputs.m_ptype);
	
	ErrOutput(u8"Affine:", u8"AffinePass Fin");
 
	if(outimgpath.size()>0)
	{
		ofstream affinefile(outimgpath + "_AffineMat.txt");
		affinefile << Affine;			
	}
				

	if(bAffineOnly)
	{
		ErrOutput(u8"Affine:","Affine Only Done quit;");
		return 0;
	}

	sWrapResults wpres;
	WrapPass(TpmInputs, tpmPriors, Affine, wpres);

	sWriteImgs sWrites;
	TpmWrite8(TpmInputs, tpmPriors, wpres, sWrites);

	if (!sWrites.FieldInverse.IsEmpty() )
	{
		sitk::WriteImage(sWrites.FieldInverse.m_Img, outimgpath + "_FieldInverse.nii");
		cout<<outimgpath + "_FieldInverse.nii"<<" Writed"<<endl;
	}
	//if (std::find(modes.begin(),modes.end(),"FieldForwared")!=modes.end())
	if(!sWrites.FieldForwared.IsEmpty())
	{
		sitk::WriteImage(sWrites.FieldForwared.m_Img, outimgpath + "_FieldForwared.nii");
		cout<<outimgpath + "_FieldForwared.nii"<<" Writed"<<endl;
	}
	
	//if (std::find(modes.begin(),modes.end(),"BiasCorrected")!=modes.end())
	//if(sWrites.)	
	for (int ichan=0;ichan<sWrites.chans.size();ichan++)
	{
		auto & chan = sWrites.chans[ichan];
		if (!chan.m_Nc.IsEmpty())
		{
			auto outfile = FormatedString("%s_BiasCorrected_%d.nii",outimgpath.c_str(),ichan);
			sitk::WriteImage(chan.m_Nc.m_Img, outfile);
			cout<<outfile<<" Writed"<<endl;
		}
	}
	{
		for (int ichan=0;ichan<sWrites.chans.size();ichan++)
		{
			auto & chan = sWrites.chans[ichan];
			if (!chan.m_Nf.IsEmpty())
			{
				auto outfile = FormatedString("%s_BiasField_%d.nii",outimgpath.c_str(),ichan);
				sitk::WriteImage(chan.m_Nf.m_Img, outfile);				
				cout<<outfile<<" Writed"<<endl;
			}
		}		
	}	
	{
		for (int ichan=0;ichan<sWrites.tissImported.size();ichan++)
		{
			auto & ti = sWrites.tissImported[ichan];
			if (!ti.IsEmpty())
			{
				auto outfile = FormatedString("%s_TissImported_%d.nii",outimgpath.c_str(),ichan);
				sitk::WriteImage(ti.m_Img, outfile);
				cout<<outfile<<" Writed"<<endl;
			}
		}		
	}	
	{
		for (int ichan=0;ichan<sWrites.tissNative.size();ichan++)
		{
			auto & tn = sWrites.tissNative[ichan];
			if (!tn.IsEmpty())
			{
				//imgs.push_back(chan.m_Img);
				sitk::WriteImage(tn.m_Img, FormatedString("%s_TissNative_%d.nii",outimgpath.c_str(),ichan));
			}
		}		
	}

	//if (std::find(modes.begin(),modes.end(),"TissModulate")!=modes.end())
	{
		for (int ichan=0;ichan<sWrites.tissModulate.size();ichan++)
		{
			auto & tm = sWrites.tissModulate[ichan];
			if (!tm.IsEmpty())
			{
				//imgs.push_back(chan.m_Img);
				sitk::WriteImage(tm.m_Img, FormatedString("%s_TissModulate_%d.nii",outimgpath.c_str(),ichan));
			}
		}		
	}
	//if (std::find(modes.begin(),modes.end(),"TissUnModulate")!=modes.end())
	{
		for (int ichan=0;ichan<sWrites.tissUnModulate.size();ichan++)
		{
			auto & tu = sWrites.tissUnModulate[ichan];
			if (!tu.IsEmpty())
			{
				//imgs.push_back(chan.m_Img);
				sitk::WriteImage(tu.m_Img, FormatedString("%s_TissUnModulate_%d.nii",outimgpath.c_str(),ichan));
			}
		}		
	}	
	return 0;
}
